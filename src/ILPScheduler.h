#ifndef CE_ILPNONPREEMPSCHEDULER_H
#define CE_ILPNONPREEMPSCHEDULER_H

#include <gurobi_c++.h>
#include <tuple>
#include <unordered_map>

#include "Instance.h"
#include "Scheduler.h"

class ILPScheduler : public Scheduler {

public:

    ILPScheduler(const Instance &inst,
                 double time_limit,
                 double target = 0,
                 double lb_time_limit = 0) :
    Scheduler(inst, time_limit),
    target(target),
    lb_time_limit(lb_time_limit) {}

    void solve() override;

private:

    double target;

    double lb_time_limit;

    bool checkSolution(std::vector< std::vector< std::vector< std::tuple <size_t, size_t,double> > > > &job_schedule);

};

#endif //CE_ILPNONPREEMPSCHEDULER_H
