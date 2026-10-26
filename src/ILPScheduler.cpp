//
// Created by bruno on 29/04/24.
//

#include "ILPScheduler.h"

#include <string>
#include <tuple>
#include <algorithm>
#include <set>
#include <chrono>
#include <gurobi_c++.h>

#include "Instance.h"
#include "WarmStart.h"
#include "util.h"

using namespace std;


class BoundCutCallback : public GRBCallback {
    GRBVar f;
    double lb;
public:
    BoundCutCallback(GRBVar f_, double lb_) : f(f_), lb(lb_) {}
protected:
    void callback() override {
        if (where != GRB_CB_MIPNODE) return;
        
        // Restrict cut generation strictly to the root node (Node 0)
        double nodecnt = getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
        if (nodecnt > 0) return; 

        if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) return;
        
        double f_lp = getNodeRel(f);         
        if (f_lp < lb - 1e-6) {
            addCut(f >= lb);                 
            cout << "Root node bound cut applied: f >= " << lb << " (LP was " << f_lp << ")\n";
        }
    }
};

void ILPScheduler::solve() {
    try {
        GRBEnv env(true);
        // env.set(GRB_IntParam_OutputFlag, 1);
        // env.set(GRB_IntParam_LogToConsole, 1);
        env.start();

        GRBModel model(env);
        model.set(GRB_DoubleParam_BestObjStop, target);

        // s_uj: Binary = 1 if sporadic task 'u' is assigned to processor 'j'
        vector<vector<GRBVar>> s_uj(inst.sporadic_tasks.size());
        vector<size_t> sporadic_prio_order(inst.sporadic_tasks.size());

        for (size_t u = 0; u < inst.sporadic_tasks.size(); ++u) {
            sporadic_prio_order[u] = u;
            s_uj[u] = vector<GRBVar>(inst.processors);
            GRBLinExpr mapping_constr_expr;

            for (size_t j = 0; j < inst.processors; ++j) {
                s_uj[u][j] =
                    model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "s_" + to_string(u) + "_" + to_string(j));
                mapping_constr_expr += s_uj[u][j];
            }

            // Constraint (3): Each sporadic task must be mapped to exactly one processor.
            model.addConstr(mapping_constr_expr == 1, "Sporadic_Map_Task_" + to_string(u));
        }

        // Sort sporadic tasks by deadline (decreasingly) for priority assignment in TDA-like constraints
        sort(sporadic_prio_order.begin(), sporadic_prio_order.end(),
             [&](size_t u, size_t v) {
                 return inst.sporadic_tasks[u].deadline - inst.sporadic_tasks[u].jitter
                            > inst.sporadic_tasks[v].deadline - inst.sporadic_tasks[v].jitter;
             });

        // This section implements constraints based on TDA
        // for each sporadic task, considering interference from higher-priority sporadic tasks.

        for (size_t u_prio_rank = 0; u_prio_rank < sporadic_prio_order.size(); ++u_prio_rank) {
            const size_t u_idx = sporadic_prio_order[u_prio_rank];
            const auto& u = inst.sporadic_tasks[u_idx];

            // Determine critical instants for the response time analysis window considering jitter
            set<double> critical_instants;

            for (size_t v_prio_rank = u_prio_rank + 1; v_prio_rank < sporadic_prio_order.size(); ++v_prio_rank) {
                const size_t v_idx = sporadic_prio_order[v_prio_rank];
                const auto& v = inst.sporadic_tasks[v_idx];

                for (int k = 1; v.min_inter_arrival * k - v.jitter <= u.deadline - u.jitter; ++k) {
                    double t = v.min_inter_arrival * k - v.jitter;
                    if (t > 0) {
                        critical_instants.insert(t);
                    }
                }
            }
            critical_instants.insert(u.deadline - u.jitter); // Always include the effective deadline in the critical instants

            // w_ut_exprs[t] holds the expression for the total demand at critical instant 't'
            vector<GRBLinExpr> w_ut_exprs(critical_instants.size(), inst.sporadic_tasks[u_idx].WCET);

            // Calculate interference from lower-priority sporadic tasks (t2)
            for (size_t v_prio_rank = u_prio_rank + 1; v_prio_rank < sporadic_prio_order.size(); ++v_prio_rank) {
                const size_t v_idx = sporadic_prio_order[v_prio_rank];
                const auto& v = inst.sporadic_tasks[v_idx];

                for (size_t j = 0; j < inst.processors; ++j) {
                    // o_uvj: Auxiliary variable indicating if both u and v are on processor j
                    GRBVar o_uvj = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                                "o_" + to_string(u_idx) + "_" + to_string(v_idx) + "_" + to_string(j));

                    // Constraints (6) and (7)
                    model.addConstr(o_uvj <= s_uj[u_idx][j], "TDA_Aux_7_" + to_string(u_idx));
                    model.addConstr(o_uvj <= s_uj[v_idx][j], "TDA_Aux_8_" + to_string(v_idx));
                    model.addConstr(o_uvj >= s_uj[u_idx][j] + s_uj[v_idx][j] - 1,
                                    "TDA_Aux_9_" + to_string(u_idx) + "_" + to_string(v_idx));

                    size_t t = 0;
                    for (double critical_time : critical_instants) {
                        // Add interference term: Demand * o_iuj
                        w_ut_exprs[t] +=
                            ceil((critical_time + v.jitter) / v.min_inter_arrival) * v.WCET * o_uvj;
                        t++;
                    }
                }
            }

            double bigM = u.WCET;

            // Iterate over the interfering tasks (v < u)
            for (size_t v_prio_rank = u_prio_rank + 1; v_prio_rank < sporadic_prio_order.size(); ++v_prio_rank) {
                const size_t v_idx = sporadic_prio_order[v_prio_rank];
                const auto& v = inst.sporadic_tasks[v_idx];

                double demand = ceil((u.deadline - u.jitter + v.jitter) / v.min_inter_arrival) * v.WCET;
                bigM += demand;
            }

            // Apply the TDA feasibility test for u
            vector<GRBVar> z_ut(critical_instants.size());
            GRBLinExpr slack_sum_expr;
            size_t instant_idx = 0;

            for (double critical_time : critical_instants) {
                z_ut[instant_idx] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                                       "z_ut_" + to_string(u_idx) + "_" + to_string(instant_idx));

                // Constraint (9): TDA feasibility check
                model.addConstr(w_ut_exprs[instant_idx] <=
                    (1 - z_ut[instant_idx]) * bigM + critical_time,
                    "TDA_Test_Task_" + to_string(u_idx) + "_Inst_" + to_string(instant_idx));

                slack_sum_expr += z_ut[instant_idx];
                instant_idx++;
            }

            // Constraint (10): At least one instant must satisfy the TDA test
            model.addConstr(slack_sum_expr == 1, "TDA_Feasibility_Task_" + to_string(u_idx));
        }

        // x_ihjk: Binary = 1 if periodic task 'i', job 'h' is assigned to processor 'j' in frame 'k'
        vector<vector<vector<vector<GRBVar>>>> x_ihjk(inst.periodic_tasks.size());

        for (size_t i_idx = 0; i_idx < inst.periodic_tasks.size(); ++i_idx) {
            const auto& i = inst.periodic_tasks[i_idx];
            x_ihjk[i_idx] = vector< vector < vector <GRBVar> > >(inst.periodic_tasks[i_idx].jobs, vector< vector <GRBVar> >(inst.processors));

            for (size_t h = 0; h < i.jobs; ++h) {
                for (size_t j = 0; j < inst.processors; ++j) {
                    x_ihjk[i_idx][h][j] = vector <GRBVar>(inst.frames);
                    for (size_t k : i.schedulable_frames[h]) {
                        x_ihjk[i_idx][h][j][k] =
                            model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                         "x_" + to_string(i_idx) + "_" + to_string(h) + "_" + to_string(j) + "_" + to_string(k));
                    }
                }
            }
        }

        // y_ij: Binary variable = 1 if periodic task 'i' is assigned to processor 'j'
        vector<vector<GRBVar>> y_ij(inst.periodic_tasks.size());
        for (size_t i_idx = 0; i_idx < inst.periodic_tasks.size(); ++i_idx) {
            y_ij[i_idx] = vector<GRBVar>(inst.processors);
            GRBLinExpr mapping_constr_expr;
            for (size_t j = 0; j < inst.processors; ++j) {
                y_ij[i_idx][j] =
                    model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + to_string(i_idx) + "_" + to_string(j));
                mapping_constr_expr += y_ij[i_idx][j];
            }
            // Constraint (2): Each periodic task must be mapped to exactly one processor.
            model.addConstr(mapping_constr_expr == 1, "Periodic_Map_Task_" + to_string(i_idx));
        }

        // Constraint (4): Lnks the partitioning variable (y_ij) with the scheduling variables (x_ihjk).
        for (size_t i_idx = 0; i_idx < inst.periodic_tasks.size(); ++i_idx) {
            const auto& i = inst.periodic_tasks[i_idx];
            for (size_t h = 0; h < i.jobs; ++h) {
                for (size_t j = 0; j < inst.processors; ++j) {
                    GRBLinExpr assignment_sum_expr = 0;
                    for (size_t f : i.schedulable_frames[h]) {
                        assignment_sum_expr += x_ihjk[i_idx][h][j][f];
                    }
                    model.addConstr(assignment_sum_expr == y_ij[i_idx][j],
                                    "Partition_Link_T" + to_string(i_idx) + "_J" + to_string(h) + "_P" + to_string(j));
                }
            }
        }

        bool all_wcets_integer = true;
        for (const auto& task : inst.periodic_tasks) {
            if (!areAlmostEqual(abs(task.WCET), round(task.WCET))) {
                all_wcets_integer = false;
                break;
            }
        }

        if (all_wcets_integer) {
            for (const auto& task : inst.sporadic_tasks) {
                if (!areAlmostEqual(abs(task.WCET), round(task.WCET))) {
                    all_wcets_integer = false;
                    break;
                }
            }
        }

        // Run warm start computation
        WarmStart ws;
        if (lb_time_limit > 0) {
            ws = computeWS(inst, lb_time_limit, all_wcets_integer);
            // Apply warm start values to partitioning variables
            for (size_t i = 0; i < inst.periodic_tasks.size(); ++i)
                for (size_t j = 0; j < inst.processors; ++j)
                    y_ij[i][j].set(GRB_DoubleAttr_Start, ws.y[i][j]);

            for (size_t i = 0; i < inst.sporadic_tasks.size(); ++i)
                for (size_t j = 0; j < inst.processors; ++j)
                    s_uj[i][j].set(GRB_DoubleAttr_Start, ws.s[i][j]);
            time_limit += lb_time_limit - ws.solution_time;

        }
        model.set(GRB_DoubleParam_TimeLimit, time_limit);

        // frame_var represents the maximum processor occupancy (frame length)
        GRBVar frame_var;
        if (all_wcets_integer)
            frame_var = model.addVar(0, inst.major_cycle, 0.0, GRB_INTEGER, "Frame_Length");
        else
            frame_var = model.addVar(0, inst.major_cycle, 0.0, GRB_CONTINUOUS, "Frame_Length");

        // Set objective: Minimize the maximum frame length
        GRBLinExpr obj = frame_var;
        model.setObjective(obj, GRB_MINIMIZE);

        // Constraint (5): Processor demand in frame 'k' on processor 'j' must be less than or equal to frame_var (f)
        for (size_t j = 0; j < inst.processors; ++j) {
            for (size_t k = 0; k < inst.frames; ++k) {
                GRBLinExpr demand_expr = 0;

                // Demand from Periodic Tasks
                for (size_t i_idx = 0; i_idx < inst.periodic_tasks.size(); ++i_idx) {
                    const auto& i = inst.periodic_tasks[i_idx];
                    for (size_t h = 0; h < i.jobs; ++h) {
                        for (size_t i_k : i.schedulable_frames[h]) {
                             if (i_k == k) {
                                 demand_expr += x_ihjk[i_idx][h][j][k] * i.WCET;
                             }
                        }
                    }
                }

                // Demand from Sporadic Tasks
                for (size_t s_idx = 0; s_idx < inst.sporadic_tasks.size(); ++s_idx) {
                    const auto& s = inst.sporadic_tasks[s_idx];
                    demand_expr += s_uj[s_idx][j] *
                        ceil(static_cast<double>(inst.minor_cycle + s.jitter) / s.min_inter_arrival) * s.WCET;
                }

                model.addConstr(demand_expr <= frame_var,
                                "Demand_P" + to_string(j) + "_F" + to_string(k));
            }
        }

        // Constraint (11): Ensure mutual exclusion
        for (size_t i_idx = 0; i_idx < inst.periodic_tasks.size(); ++i_idx) {
            const auto& i = inst.periodic_tasks[i_idx];
            for (size_t h = 0; h < i.jobs; ++h) {

                // Get unique list of tasks that conflict with i/h:
                vector<size_t> related_tasks = i.excl_relation;
                // Remove redundant checks (only consider i < i_related_task)
                related_tasks.erase(remove_if(related_tasks.begin(), related_tasks.end(),
                                               [&](size_t i_rel_idx) { return i_idx > i_rel_idx; }), related_tasks.end());
                // Add precedence-related tasks (a task must finish before its successor starts)
                related_tasks.insert(related_tasks.end(), i.prec_relation.begin(), i.prec_relation.end());

                for (size_t i_rel_idx : related_tasks) {
                    const auto& i_rel = inst.periodic_tasks[i_rel_idx];
                    for (size_t l = 0; l < i_rel.jobs; ++l) {

                        // Find mutual frames where both jobs are schedulable
                        const vector<size_t> &k_h = i.schedulable_frames[h];
                        const vector<size_t> &k_l = i_rel.schedulable_frames[l];
                        vector<size_t> mutual_frames;
                        set_intersection(k_h.begin(), k_h.end(), k_l.begin(), k_l.end(), back_inserter(mutual_frames));

                        for (size_t k : mutual_frames) {
                            GRBLinExpr sum_x1 = 0;
                            GRBLinExpr sum_x2 = 0;

                            // Sum of assignment variables over all processors for the mutual frame
                            for (size_t proc_idx = 0; proc_idx < inst.processors; ++proc_idx) {
                                // Variable retrieval assuming correct sparse indexing
                                sum_x1 += x_ihjk[i_idx][h][proc_idx][k];
                                sum_x2 += x_ihjk[i_rel_idx][l][proc_idx][k];
                            }

                            // Constraint (11): Only one of the conflicting jobs can be scheduled in the mutual frame
                            model.addConstr(sum_x1 + sum_x2 <= 1,
                                            "Excl_T" + to_string(i_idx) + "J" + to_string(h) +
                                            "_T" + to_string(i_rel_idx) + "J" + to_string(l) + "_F" + to_string(k));
                        }
                    }
                }
            }
        }

        // Constraint (12): Precedence relations
        for (size_t i_idx = 0; i_idx < inst.periodic_tasks.size(); ++i_idx) {
            const auto& i = inst.periodic_tasks[i_idx];
            for (size_t h = 0; h < i.jobs; ++h) {

                for (size_t prec_task_idx : i.prec_relation) {
                    const vector<size_t> &k_h = i.schedulable_frames[h]; // Successor's schedulable frames
                    GRBLinExpr predecessor_sum = 0; // Sum of x for predecessor task (prec_task)
                    for (size_t k : k_h) {
                        GRBLinExpr successor_sum = 0; // Sum of x for successor task (i)
                        for (size_t j = 0; j < inst.processors; ++j) {
                            successor_sum += x_ihjk[i_idx][h][j][k];
                        }

                        model.addConstr(successor_sum <= predecessor_sum,
                                        "Preced_Succ" + to_string(i_idx) + "_Pred" + to_string(prec_task_idx) + "_J" + to_string(h) + "_F" + to_string(k));

                        for (size_t j = 0; j < inst.processors; ++j) {
                            predecessor_sum += x_ihjk[prec_task_idx][h][j][k];
                        }
                    }
                }
            }
        }        

        // Register the pre-computed bound cut callback if warm start produced a bound
        BoundCutCallback* cb = nullptr;
        if (lb_time_limit > 0 && ws.Cvalue_bound > 0) {
            cb = new BoundCutCallback(frame_var, ws.Cvalue_bound);
            model.setCallback(cb);
            cout << "BoundCutCallback registered with pre-computed bound: " << ws.Cvalue_bound << "\n";
        }

        cout << "Starting ILP computation with "
                  << inst.periodic_tasks.size() << " periodic tasks and "
                  << inst.sporadic_tasks.size() << " sporadic tasks on "
                  << inst.processors << " processors.\n";
        model.optimize();
        cout << "ILP computation completed.\n";
        double gap = model.get(GRB_DoubleAttr_MIPGap);

        solution_time = model.get(GRB_DoubleAttr_Runtime) + ws.solution_time;

        if (model.get(GRB_IntAttr_SolCount) > 0) {
            max_occupancy = model.get(GRB_DoubleAttr_ObjVal);
            cout << "Maximum Frame Occupancy (f) = " << max_occupancy << " GAP = " << gap  << endl;
            cout << "Total solution time (s) = " << solution_time << endl;

            // Extract the schedule and assignment from the solution variables

            // Reinitialize structures for the solution
            // job_schedule[task][job] -> list of (processor, frame, WCET)
            vector<vector<vector<tuple<size_t, size_t, double>>>> job_schedule(inst.periodic_tasks.size());
            for (size_t i = 0; i < inst.periodic_tasks.size(); ++i) {
                job_schedule[i].resize(inst.periodic_tasks[i].jobs);
            }

            // Extract periodic task schedule
            for (size_t frame = 0; frame < inst.frames; ++frame) {
                for (size_t task = 0; task < inst.periodic_tasks.size(); ++task) {
                    for (size_t job = 0; job < inst.periodic_tasks[task].jobs; ++job) {
                        for (size_t processor = 0; processor < inst.processors; ++processor) {
                            // Check variable existence and solution value
                            if (x_ihjk[task][job][processor][frame].index() >= 0 &&
                                x_ihjk[task][job][processor][frame].get(GRB_DoubleAttr_X) > 0.5) {
                                frame_schedule[frame][processor].emplace_back(task, job, inst.periodic_tasks[task].WCET);
                                job_schedule[task][job].emplace_back(processor, frame, inst.periodic_tasks[task].WCET);
                            }
                        }
                    }
                }
            }

            // Extract sporadic task assignment
            for (size_t task = 0; task < inst.sporadic_tasks.size(); ++task) {
                for (size_t processor = 0; processor < inst.processors; ++processor) {
                    if (s_uj[task][processor].index() >= 0 &&
                        s_uj[task][processor].get(GRB_DoubleAttr_X) > 0.5) {
                        // Store: (processor, priority_rank, WCET)
                        sporadic_task_assignment[inst.sporadic_tasks[task].id] =
                            make_tuple(processor, task, inst.sporadic_tasks[task].WCET);
                        break; // Stop after finding the assigned processor
                    }
                }
            }

            valid_solution = checkSolution(job_schedule);
            if (!valid_solution) {
                cout << "Warning: The extracted solution does not satisfy all constraints." << endl;
                exit(1);
            }


            cout << "\n======================================================\n";
            cout << " CYCLIC EXECUTIVE SCHEDULE (FRAME-BY-FRAME)  \n";
            cout << "======================================================\n";

            for (size_t frame = 0; frame < frame_schedule.size(); ++frame) {
                cout << "\n------------------------------------------------------\n";
                cout << "**FRAME " << frame << "**: Minor Cycle Time Slot\n";
                cout << "------------------------------------------------------\n";

                for (size_t processor = 0; processor < inst.processors; ++processor) {
                    cout << "  [Processor " << processor << "]: ";

                    const auto& assignments = frame_schedule[frame][processor];

                    if (assignments.empty()) {
                        cout << "IDLE\n";
                    } else {
                        // Iterate through assignments for this (processor, frame) slot
                        for (const auto& assignment : assignments) {
                            size_t task_id = get<0>(assignment);
                            size_t job_id = get<1>(assignment);
                            double wcet = get<2>(assignment);

                            cout << "Periodic Task " << task_id << " (Job " << job_id
                                      << ") [WCET: " << wcet << "]\n";
                        }
                    }
                }
            }

            cout << "\n======================================================\n";
            cout << " SPORADIC TASK PROCESSOR ASSIGNMENT \n";
            cout << "======================================================\n";

            cout << "Task ID    | Core      |  Priority | WCET\n";
            cout << "-----------|-----------|-----------|------\n";

            for (size_t task = 0; task < inst.sporadic_tasks.size(); ++task) {
                size_t task_id = inst.sporadic_tasks[task].id;
                // Check if the assignment map contains the task_id
                if (task_id < sporadic_task_assignment.size()) {
                    const auto& assignment = sporadic_task_assignment.at(task_id);
                    size_t assigned_processor = get<0>(assignment);
                    size_t priority = get<1>(assignment);
                    double wcet = get<2>(assignment); // Use double to match the tuple type

                    cout << "S_Task " << task_id << "   | "
                              << assigned_processor << "        | "
                              << priority << "         | "
                              << wcet << "\n";
                }
            }

            cout << "\n======================================================\n";
        } else {
            cout << "No feasible solution found within the specified time limit or constraints." << endl;
            valid_solution = false;
        }

    } catch (GRBException& e) {
        cerr << "Gurobi Error: " << e.getMessage() << " (Code: " << e.getErrorCode() << ")" << endl;
    } catch (exception& e) {
        cerr << "Standard Exception: " << e.what() << endl;
    }
}
bool ILPScheduler::checkSolution(vector< vector< vector< tuple <size_t, size_t,double> > > > &job_schedule) {
    vector< vector<double> > frame_occupancy(inst.processors, vector<double>(inst.frames, 0));
    for (size_t f = 0; f < inst.frames; f++) {
        for (size_t p = 0; p < inst.processors; p++) {
            for (size_t s = 0; s < frame_schedule[f][p].size(); s++) {
                size_t t = get<0>(frame_schedule[f][p][s]);
                size_t j = get<1>(frame_schedule[f][p][s]);
                double time = get<2>(frame_schedule[f][p][s]);
                frame_occupancy[p][f] += time;

                if (not(inst.periodic_tasks[t].schedulable_frames[j].front() <= f &&
                        f <= inst.periodic_tasks[t].schedulable_frames[j].back())) {
                    return false;
                }

                for (size_t t2: inst.periodic_tasks[t].prec_relation) {
                    if (get<1>(job_schedule[t][j].back()) < get<1>(job_schedule[t2][j].front()) ||
                        (get<1>(job_schedule[t][j].back()) == get<1>(job_schedule[t2][j].front())
                         && get<0>(job_schedule[t][j].back()) != get<0>(job_schedule[t2][j].front())))
                        return false;
                }

                for (size_t t2: inst.periodic_tasks[t].succ_relation) {
                    if (get<1>(job_schedule[t][j].back()) > get<1>(job_schedule[t2][j].front()) ||
                        (get<1>(job_schedule[t][j].back()) == get<1>(job_schedule[t2][j].front())
                         && get<0>(job_schedule[t][j].back()) != get<0>(job_schedule[t2][j].front())))
                        return false;
                }

                vector<size_t> j1_frames;
                for (auto& schedule : job_schedule[t][j]) {
                    j1_frames.push_back(std::get<1>(schedule));
                }
                for (size_t t2: inst.periodic_tasks[t].excl_relation) {
                    for (size_t j2 = 0; j2 < inst.periodic_tasks[t2].jobs; j2++) {
                        vector<size_t> j2_frames;
                        for (auto& schedule : job_schedule[t2][j2]) {
                            j2_frames.push_back(std::get<1>(schedule));
                        }
                        vector<size_t> mutual_frames;
                        set_intersection(j1_frames.begin(), j1_frames.end(), j2_frames.begin(), j2_frames.end(), back_inserter(mutual_frames));
                        for (size_t mutual_f : mutual_frames) {
                            auto it_j1 = find_if(job_schedule[t][j].begin(), job_schedule[t][j].end(),
                                                   [mutual_f](const auto& tuple) { return get<1>(tuple) == mutual_f; });
                            auto it_j2 = find_if(job_schedule[t2][j2].begin(), job_schedule[t2][j2].end(),
                                                 [mutual_f](const auto& tuple) { return get<1>(tuple) == mutual_f; });
                            if (get<0>(*it_j1) != get<0>(*it_j2)) {
                                return false;
                            }
                        }
                    }
                }
            }
        }
    }

    for (size_t p = 0; p < inst.processors; p++) {
        for (size_t f = 0; f < inst.frames; f++) {
            if (isLessThan(max_occupancy, frame_occupancy[p][f]) ) {
                max_occupancy = frame_occupancy[p][f];
            }
        }
    }
    return true;
}