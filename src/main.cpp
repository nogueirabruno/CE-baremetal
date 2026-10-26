#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <optional>
#include <boost/integer/common_factor.hpp>

#include "Instance.h"
#include "ILPScheduler.h"
#include "util.h"

using namespace std;

pair<int, int> calculateLCMGCD(const vector<int>& numbers) {
    if (numbers.empty()) {
        return std::make_pair(0, 0);
    }

    int lcm = numbers[0];
    int gcd = numbers[0];

    // Iterate through the numbers in the vector and update the lcm and gcd
    for (size_t i = 1; i < numbers.size(); ++i) {
        lcm = boost::integer::lcm(lcm, numbers[i]);
        gcd = boost::integer::gcd(gcd, numbers[i]);
    }

    return make_pair(lcm, gcd);
}

Instance read_instance(const string &file_name){
    Instance i;

    ifstream inputFile(file_name);

    if (!inputFile.is_open()) {
        throw runtime_error("Error: Unable to open the input file '" + file_name + "'");
    }

    int lineNumber = 0;
    string line;

    do {
        getline(inputFile, line);
        lineNumber++;
    } while (line[0] == '#' || line.empty());

    stringstream ss(line);
    int periodic_tasks, sporadic_tasks, processors;
    optional<size_t> frame_size;

    if (!(ss >> periodic_tasks >> sporadic_tasks >> processors)) {
        throw std::runtime_error("Error at line " + to_string(lineNumber) + ": Unable to read the number of tasks and processors in the input file '" + file_name +"'");
    }

    if (periodic_tasks < 0 || sporadic_tasks < 0 || processors <= 0) {
        throw std::runtime_error("Error: Number of tasks and processors cannot be negative.");
    }

    i.processors = processors;
    i.periodic_tasks = vector<PeriodicTask>(periodic_tasks);
    i.sporadic_tasks = vector<SporadicTask>(sporadic_tasks);
    i.total_job_relations = 0;

    int temp;
    if (ss >> temp) {
        frame_size = temp;
    }

    double max_wcet = 0;

    vector<int> periods;
    for (size_t t = 0; t < i.periodic_tasks.size() + i.sporadic_tasks.size(); t++) {
        lineNumber++;
        getline(inputFile, line);
        if (line[0] == '#' || line.empty()) {t--; continue;};

        ss.str("");
        ss.clear();
        ss << line;
        if (t < i.periodic_tasks.size()) {
            i.periodic_tasks[t].id = t;
            if (!(ss >> i.periodic_tasks[t].period >> i.periodic_tasks[t].WCET)) {
                throw std::runtime_error("Error at line " + to_string(lineNumber) + ": Unable to read the period and WCET of periodic task " +
                                         to_string(t + 1) + " in the input file '" + file_name + "'");
            }
            if (i.periodic_tasks[t].WCET > max_wcet) {
                max_wcet = i.periodic_tasks[t].WCET;
            }
            i.periodic_tasks[t].excl_relation2 = vector<bool>(periodic_tasks, false);
            i.periodic_tasks[t].prec_relation2 = vector<bool>(periodic_tasks, false);
            periods.push_back(i.periodic_tasks[t].period);
        } else {
            if (!(ss >> i.sporadic_tasks[t-periodic_tasks].min_inter_arrival >> i.sporadic_tasks[t-periodic_tasks].WCET
                >> i.sporadic_tasks[t-periodic_tasks].deadline >> i.sporadic_tasks[t-periodic_tasks].jitter)) {
                throw std::runtime_error("Error at line " + to_string(lineNumber) + ": Unable to read the minimum interarrival time and WCET of sporadic task " +
                                         to_string(t-periodic_tasks + 1) + " in the input file '" + file_name + "'");
            }
            if (i.sporadic_tasks[t-periodic_tasks].WCET > max_wcet) {
                max_wcet = i.sporadic_tasks[t-periodic_tasks].WCET;
            }
        }
    }

    tie(i.major_cycle, i.minor_cycle) = calculateLCMGCD(periods);
    if (frame_size.has_value()) {
        i.minor_cycle = *frame_size;

        if (i.major_cycle % i.minor_cycle != 0) {
            throw std::runtime_error(
                "Invalid minor cycle: provided value (" + std::to_string(i.minor_cycle) +
                ") does not divide the major cycle (" + std::to_string(i.major_cycle) + ")."
            );
        }

        if (isLessThan(i.minor_cycle, max_wcet)) {
            throw std::runtime_error(
                "Invalid minor cycle: provided value (" + std::to_string(i.minor_cycle) +
                ") is smaller than the maximum WCET among all tasks (" +
                std::to_string(max_wcet) + ")."
            );
        }
    } else {
        if (isLessThan(i.minor_cycle, max_wcet)) {
            cout <<
                "WARNING: Invalid minor cycle: computed value (" + std::to_string(i.minor_cycle) +
                ") is smaller than the maximum WCET among all tasks (" +
                std::to_string(max_wcet) + ").\n";
        }
    }
    i.frames = i.major_cycle/i.minor_cycle;

    i.total_jobs = 0;
    for (int t = 0; t < periodic_tasks; t++) {
        i.periodic_tasks[t].jobs = i.major_cycle / i.periodic_tasks[t].period;
        i.total_jobs += i.periodic_tasks[t].jobs;
        i.periodic_tasks[t].schedulable_frames = std::vector< std::vector<size_t> >(i.periodic_tasks[t].jobs);
        for (size_t j = 0; j < i.major_cycle / i.periodic_tasks[t].period; j++) {
            size_t release_time = i.periodic_tasks[t].period * j;
            size_t deadline     = i.periodic_tasks[t].period * (j + 1);

            size_t first_frame = (release_time + i.minor_cycle - 1) / i.minor_cycle;
            size_t last_frame  = ((deadline - i.minor_cycle) / i.minor_cycle) + 1;

            std::vector<size_t> schedulable_frames;
            schedulable_frames.resize(last_frame - first_frame);
            std::iota(schedulable_frames.begin(), schedulable_frames.end(), first_frame);
            i.periodic_tasks[t].schedulable_frames[j] = schedulable_frames;
        }
    }

    // Sort by min_inter_arrival (ascending)
    std::sort(i.sporadic_tasks.begin(), i.sporadic_tasks.end(),
              [](const auto& a, const auto& b) {
                  return a.deadline < b.deadline;
              });
    for(size_t t = 0; t < i.sporadic_tasks.size(); t++) {
        i.sporadic_tasks[t].id = t;
    }

    while (!inputFile.eof()) {
        lineNumber++;
        getline(inputFile, line);
        if (line[0] == '#' || line.empty()) continue;
        ss.str("");
        ss.clear();
        ss << line;

        char relation;
        int task1, task2;
        if (!(ss >> relation >> task1 >> task2)) {
            throw std::runtime_error("Error at line " + to_string(lineNumber) + ": Unable to read the task relation in the input file '" + file_name + "'");
        }
        if (task1 < 0 || task1 >= periodic_tasks || task2 < 0 || task2 >= periodic_tasks) {
            throw std::runtime_error("Error at line " + to_string(lineNumber) + ": To define an exclusion or precedence relation, tasks " + to_string(task1) + " and " + to_string(task2) + " must be periodic tasks in input file '" + file_name + "'");
        }
        if (relation == 'x') {
                i.periodic_tasks[task1].excl_relation.push_back(task2);
                i.periodic_tasks[task1].excl_relation2[task2] = true;
                i.periodic_tasks[task2].excl_relation.push_back(task1);
                i.periodic_tasks[task2].excl_relation2[task1] = true;
                i.total_job_relations += i.periodic_tasks[task1].jobs;
        } else if (relation == 'p') {
            if (i.periodic_tasks[task1].period != i.periodic_tasks[task2].period) {
                throw std::runtime_error("Error at line " + to_string(lineNumber) + ": To define a precedence relation, tasks " + to_string(task1) + " and " + to_string(task2) + " need to have the same period in input file '" + file_name + "'");
            }
            i.periodic_tasks[task2].prec_relation.push_back(task1);
            i.periodic_tasks[task2].prec_relation2[task1] = true;
            i.periodic_tasks[task1].succ_relation.push_back(task2);
            i.total_job_relations += i.periodic_tasks[task1].jobs;
        } else {
            throw std::runtime_error("Error at line " + to_string(lineNumber) + ": Invalid task relation in the input file '" + file_name + "'");
        }

    }
    inputFile.close();

    // Compute utilization
    i.utilization = 0.0;
    for (auto &task : i.periodic_tasks)
        i.utilization += task.WCET / task.period;
    for (auto &task : i.sporadic_tasks)
        i.utilization += task.WCET / task.min_inter_arrival;

    return i;
}

int main(int argc, char* argv[]) {

    try {
        string file_name;
        int time_limit;
        double target = 0.0;
        int ws_time_limit = 0;   // Only defined if -lb is used

        if (argc < 3) {
            throw invalid_argument(
                "Usage: " + string(argv[0]) +
                " <file_name> <time_limit> [-ws <ws_time>] [-T <target>]"
            );
        }

        // Mandatory arguments
        file_name = argv[1];

        try {
            time_limit = stoi(argv[2]);
        } catch (...) {
            throw invalid_argument("Invalid time limit: must be an integer.");
        }

        // Parse optional flags
        for (int i = 3; i < argc; ++i) {
            string arg = argv[i];

            if (arg == "-ws") {
                // Check next argument exists
                if (i + 1 >= argc)
                    throw invalid_argument("Missing <lb_time> after -ws.");

                // Try to parse time
                try {
                    ws_time_limit = stod(argv[++i]);
                } catch (...) {
                    throw invalid_argument("Invalid numeric value for -ws.");
                }
            }
            else if (arg == "-T") {
                if (++i >= argc)
                    throw invalid_argument("Missing <target> value for -T.");

                try {
                    target = stod(argv[i]);
                } catch (...) {
                    throw invalid_argument("Invalid numeric value for -T.");
                }
            }
            else {
                throw invalid_argument("Unknown option: " + arg);
            }
        }

        cout << "\nReading instance " << file_name << "...";
        Instance inst = read_instance(file_name);
        cout << " Done!\n";

        // Instance information
        cout << "Instance information:\n";
        cout << "  Periodic tasks: " << inst.periodic_tasks.size()
             << ", Jobs: " << inst.total_jobs
             << ", Job relations: " << inst.total_job_relations << "\n";
        cout << "  Sporadic tasks: " << inst.sporadic_tasks.size() << "\n";
        cout << "  Processors: " << inst.processors
             << ", Frames: " << inst.frames
             << ", Minor cycle: " << inst.minor_cycle
             << ", Major cycle: " << inst.major_cycle << "\n";
        cout << "  System utilization: " << inst.utilization << "\n\n";

        if (ws_time_limit > 0) {
            cout << "Warm start ENABLED with time limit = "
                 << ws_time_limit << "s\n\n";
        }

        if (target > 0) {
            cout << "Target solution ENABLED with value = "
                 << target << "s\n\n";
        }

        cout << "Starting ILP Scheduler with "
             << time_limit << "s main time limit...\n\n";

        Scheduler* sched = new ILPScheduler(inst, time_limit, target, ws_time_limit);
        sched->solve();
        delete sched;

        return 0;

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}