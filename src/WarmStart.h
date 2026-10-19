//
// Created by optlab on 23/11/25.
//

#ifndef WARMSTART_H
#define WARMSTART_H

#include <gurobi_c++.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

#include "Instance.h"
#include "util.h"

struct WarmStart {
    double Cvalue;
    double Cvalue_bound;
    bool optimal;
    double solution_time;
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> s;
};

class PackingConstraintCallback : public GRBCallback {

public:

    const Instance& inst;

    std::vector<std::vector<GRBVar>>& y;
    std::vector<std::vector<GRBVar>>& s;

    GRBVar& frame_var;

    WarmStart &ws;

    PackingConstraintCallback(const Instance& inst,
                     std::vector<std::vector<GRBVar>>& y,
                     std::vector<std::vector<GRBVar>>& s,
                     GRBVar& C,
                     WarmStart &ws)
        : inst(inst), y(y), s(s), frame_var(C), ws(ws) {}

protected:

    void callback() override {

        if (where != GRB_CB_MIPSOL) return;

        double Cval = getSolution(frame_var);
        size_t J = inst.processors;

        // Define threshold divisors
        std::vector<double> thresholds = {2,3};

        // Iterate through all processors to FIND a violation
        for (size_t j = 0; j < J; j++) {
            std::vector<int> periodic_on_j;
            std::vector<int> periodic_on_j_violating;
            std::vector<int> sporadic_on_j;

            for (size_t i = 0; i < inst.periodic_tasks.size(); i++) {
                if (getSolution(y[i][j]) > 0.5) {
                    periodic_on_j.push_back(i);
                }
            }

            for (size_t u = 0; u < inst.sporadic_tasks.size(); u++) {
                if (getSolution(s[u][j]) > 0.5) {
                    sporadic_on_j.push_back(u);
                }
            }

            // Compute sporadic sum for this processor
            double sporadic_sum = 0.0;
            for (int i : sporadic_on_j) {
                sporadic_sum += inst.sporadic_tasks[i].WCET;
            }

            // Test each threshold
            for (double threshold : thresholds) {
                std::vector<double> large_wcets;
                int job_count = 0;
                periodic_on_j_violating.clear();
                // Collect jobs with effective WCET > C/threshold
                for (int i : periodic_on_j) {
                    double effective_wcet = inst.periodic_tasks[i].WCET + sporadic_sum;
                    double threshold_value = Cval / threshold;

                    if (isLessThan(threshold_value, effective_wcet)) {
                        job_count += inst.periodic_tasks[i].jobs;
                        periodic_on_j_violating.push_back(i);
                        for (size_t k = 0; k < inst.periodic_tasks[i].jobs; k++) {
                            large_wcets.push_back(effective_wcet);
                        }
                    }
                }

                // Check for violation: if more than (threshold-1)*frames jobs are large
                // E.g., for threshold=2 (C/2): violation if more than F jobs
                //       for threshold=3 (C/3): violation if more than 2F jobs
                int max_allowed_jobs = (int)std::floor((threshold - 1.0) * inst.frames);

                if (job_count > max_allowed_jobs) {
                    // Violation found: compute C_min
                    std::sort(large_wcets.begin(), large_wcets.end(), std::greater<double>());

                    double W_critical = large_wcets[max_allowed_jobs];
                    double C_min_required = threshold * W_critical;

                    // Build cut
                    GRBLinExpr sum_assignment_change = 0;

                    for (int i : periodic_on_j_violating) {
                        sum_assignment_change += (1.0 - y[i][j]);
                    }

                    for (int i : sporadic_on_j) {
                        sum_assignment_change += (1.0 - s[i][j]);
                    }

                    // Add lazy constraint
                    GRBLinExpr cut_rhs = C_min_required - C_min_required * sum_assignment_change;
                    addLazy(frame_var >= cut_rhs);
                    return;
                }
            }
        }
    }
};

WarmStart computeWS(const Instance &inst, double lb_time_limit, bool all_integers)
{
    GRBEnv env(true);
    env.set(GRB_IntParam_OutputFlag, 0);
    env.set(GRB_IntParam_LogToConsole, 0);
    env.start();

    GRBModel model(env);
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_TimeLimit, lb_time_limit);


    std::vector<std::vector<GRBVar>> y_ij(inst.periodic_tasks.size(),
                                        std::vector<GRBVar>(inst.processors));
    std::vector<std::vector<GRBVar>> s_uj(inst.sporadic_tasks.size(),
                                        std::vector<GRBVar>(inst.processors));

    for (size_t i = 0; i < inst.periodic_tasks.size(); i++)
        for (size_t j = 0; j < inst.processors; j++)
            y_ij[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

    for (size_t u = 0; u < inst.sporadic_tasks.size(); u++)
        for (size_t j = 0; j < inst.processors; j++)
            s_uj[u][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

    // frame lenght variable
    GRBVar frame_var;
    if (all_integers)
        frame_var = model.addVar(0.0, inst.major_cycle, 1.0, GRB_INTEGER, "C");
    else
        frame_var = model.addVar(0.0, inst.major_cycle, 1.0, GRB_CONTINUOUS, "C");

    // Constraint (20)
    for (size_t i = 0; i < inst.periodic_tasks.size(); i++) {
        GRBLinExpr sum = 0;
        for (size_t j = 0; j < inst.processors; j++)
            sum += y_ij[i][j];
        model.addConstr(sum == 1.0);
    }

    // Constraint (21)
    for (size_t u = 0; u < inst.sporadic_tasks.size(); u++) {
        GRBLinExpr sum = 0;
        for (size_t j = 0; j < inst.processors; j++)
            sum += s_uj[u][j];
        model.addConstr(sum == 1.0);
    }

    // Constraint (22)
    for (size_t j = 0; j < inst.processors; j++) {
        GRBLinExpr load = 0;
        for (size_t i = 0; i < inst.periodic_tasks.size(); i++)
            load += inst.periodic_tasks[i].WCET/inst.periodic_tasks[i].period * inst.minor_cycle * y_ij[i][j];
        for (size_t u = 0; u < inst.sporadic_tasks.size(); u++)
            load += inst.sporadic_tasks[u].WCET *
                ceil((inst.minor_cycle + inst.sporadic_tasks[u].jitter) / inst.sporadic_tasks[u].min_inter_arrival) * s_uj[u][j];
        model.addConstr(load <= frame_var);
    }

    // Constraint (23)
    for (size_t j = 0; j < inst.processors; j++) {
        GRBLinExpr load = 0;
        for (size_t u = 0; u < inst.sporadic_tasks.size(); u++)
            load += inst.sporadic_tasks[u].WCET *
                ceil((inst.minor_cycle + inst.sporadic_tasks[u].jitter) / inst.sporadic_tasks[u].min_inter_arrival) * s_uj[u][j];
        for (size_t i = 0; i < inst.periodic_tasks.size(); i++) {
            model.addConstr(load + inst.periodic_tasks[i].WCET * y_ij[i][j] <= frame_var);
        }
    }

    // Minimize frame lenght
    std::cout << "Starting WS computation with "
              << inst.periodic_tasks.size() << " periodic tasks and "
              << inst.sporadic_tasks.size() << " sporadic tasks on "
              << inst.processors << " processors.\n";

    WarmStart ws;
    PackingConstraintCallback cb_jl(inst, y_ij, s_uj, frame_var, ws);
    model.setCallback(&cb_jl);
    model.set(GRB_IntParam_LazyConstraints, 1);

    GRBLinExpr obj = frame_var;
    model.setObjective(obj, GRB_MINIMIZE);
    model.optimize();

    ws.optimal = (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL);
    ws.Cvalue_bound = model.get(GRB_DoubleAttr_ObjBound);
    ws.Cvalue = model.get(GRB_DoubleAttr_ObjVal);

    ws.Cvalue = frame_var.get(GRB_DoubleAttr_X);

    ws.y.resize(inst.periodic_tasks.size(),
                std::vector<double>(inst.processors, 0.0));
    ws.s.resize(inst.sporadic_tasks.size(),
                std::vector<double>(inst.processors, 0.0));

    for (size_t i = 0; i < inst.periodic_tasks.size(); i++)
        for (size_t j = 0; j < inst.processors; j++)
            ws.y[i][j] = y_ij[i][j].get(GRB_DoubleAttr_X);

    for (size_t i = 0; i < inst.sporadic_tasks.size(); i++)
        for (size_t j = 0; j < inst.processors; j++)
            ws.s[i][j] = s_uj[i][j].get(GRB_DoubleAttr_X);

    ws.solution_time = model.get(GRB_DoubleAttr_Runtime);

    std::cout << "WS computation completed in " << ws.solution_time << " seconds.\n";
    std::cout << "WS value: " << ws.Cvalue << ", Bound: " << ws.Cvalue_bound << ", Optimal: " << (ws.optimal ? "Yes" : "No") << "\n\n";

    return ws;
}

#endif //WARMSTART_H