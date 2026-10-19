#ifndef INSTANCE_H
#define INSTANCE_H

#include <vector>
#include <cmath>
#include <numeric>

struct PeriodicTask {
	size_t id;
	int period;
	double WCET;
    size_t jobs;
    std::vector< std::vector<size_t> > schedulable_frames;
    std::vector<size_t> excl_relation;
	std::vector<size_t> prec_relation;
    std::vector<size_t> succ_relation;
    std::vector<bool> excl_relation2;
    std::vector<bool> prec_relation2;
};

struct SporadicTask {
    size_t id;
    double min_inter_arrival;
    double deadline;
    double WCET;
    int priority;
    double jitter;
    std::vector<double> critical_instants;
};

struct Instance {

    size_t major_cycle, minor_cycle;

    size_t processors;

    size_t frames;

    std::vector<PeriodicTask> periodic_tasks;

    std::vector<SporadicTask> sporadic_tasks;

    size_t total_jobs;

    size_t total_job_relations;

    double utilization;

    bool isJobFrame(size_t t, size_t j, size_t f) {
        size_t first_frame = floor(((float) periodic_tasks[t].period * j)/minor_cycle);
        size_t last_frame =  floor(((float) periodic_tasks[t].period * (j + 1))/minor_cycle);

        return (first_frame <= f && f <= last_frame);
    }
    
};

#endif
