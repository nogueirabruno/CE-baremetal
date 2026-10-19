#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <vector>
#include "Instance.h"

class Scheduler {

public:

    Instance const &inst;

    double time_limit;

    bool valid_solution;

    std::vector< std::vector< std::vector< std::tuple<size_t, size_t, double> >>> frame_schedule;

    std::vector<std::tuple<size_t,size_t,double>> sporadic_task_assignment;

    double solution_time;

    double max_occupancy;

    Scheduler(const Instance &inst, double time_limit) : inst(inst), time_limit(time_limit),
            valid_solution(false), frame_schedule(inst.frames, std::vector< std::vector<
            std::tuple<size_t, size_t, double> >>(inst.processors)),
            sporadic_task_assignment(inst.sporadic_tasks.size()) {
    }

    virtual ~Scheduler() = default;

    virtual void solve() = 0;

};

#endif //SCHEDULER_H
