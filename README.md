# CE-baremetal

Official implementation and benchmark instances for the paper:

> **Unified Scheduling of Periodic and Sporadic Tasks in Multicore Bare-Metal Systems via MILP**  
> Bruno Nogueira, Pedro Rocha, Cleiber Silva Junnior, Alfredo Lima, Rodrigo Paes, and Francisco Airton Silva  
> *Journal of Systems Architecture (JSA)*

## Overview

This repository provides a Mixed-Integer Linear Programming (MILP) framework for the automated synthesis of **Cyclic Executive (CE)** schedules for multicore bare-metal real-time embedded systems. The framework jointly partitions and schedules:
- **Periodic tasks** executed in the main super loop, subject to mutual-exclusion and equal-period precedence constraints.
- **Interrupt-driven sporadic tasks** (ISRs) whose schedulability is certified offline through an embedded Time-Demand Analysis (TDA) with release jitter.

To accelerate convergence on larger instances, the solver also includes a **relaxed Bin Packing-based formulation (`R-MILP`)** with dynamic lazy constraints that computes tight lower bounds and generates warm-start task-to-core assignments for the exact MILP (`R-MILP -> MILP`).

## Project Structure

```text
CE-baremetal/
├── Makefile             # Build configuration (C++17 + Gurobi)
├── src/
│   ├── main.cpp         # Entry point, CLI parser, and instance reader
│   ├── Instance.h       # Data structures for periodic/sporadic tasks & instances
│   ├── Scheduler.h      # Base scheduler interface
│   ├── ILPScheduler.h   # Gurobi MILP scheduler declaration
│   ├── ILPScheduler.cpp # Exact MILP formulation, TDA constraints, and solution checker
│   ├── WarmStart.h      # Relaxed BPP model (R-MILP) for lower bounds and warm start
│   └── util.h           # Floating-point and helper utilities
├── instances/           # Benchmark instances evaluated in the paper
│   ├── generic/         # Generic benchmark set (generic_30, generic_50, generic_100)
│   ├── vehicle/         # Autonomous vehicle benchmark set (vehicle1, vehicle2, vehicle3)
│   ├── automotive/      # Automotive benchmark set (inst1, inst2, inst3)
│   └── t*/              # UUniFast synthetic benchmark sets
└── scalability_data/    # Sensitivity/scalability instances varying frames (f), cores (p), and utilization (u)
```

## Prerequisites

- **Gurobi Optimizer** (C++ API and shared libraries)

## Building

Run `make` from the repository root, overriding `GUROBI_HOME` and `GUROBI_VER` to match your local Gurobi installation:

```bash
make GUROBI_HOME=/path/to/gurobi/linux64 GUROBI_VER=120
```

This compiles the source files in `src/` and produces the `ce_scheduler` binary in the root directory.

To clean build artifacts:

```bash
make clean
```

## Usage

```bash
./ce_scheduler <file_name> <time_limit> [-lb <lb_time>] [-T <target>]
```

### Arguments
- `<file_name>`: Path to the input instance file (e.g., `instances/automotive/inst1`).
- `<time_limit>`: Maximum total solver runtime in seconds (integer).
- `-lb <lb_time>` *(optional)*: Enable the two-stage relaxation-guided strategy (`R-MILP -> MILP`), allocating up to `<lb_time>` seconds to the relaxed model (`R-MILP`) to compute a lower bound and warm-start assignment before solving the full MILP.
- `-T <target>` *(optional)*: Stop the solver early once an objective value of `<target>` is reached (`BestObjStop`).

### Examples

Solve the full `MILP` directly with a 300-second time limit:
```bash
./ce_scheduler instances/automotive/inst1 300
```

Solve using the two-stage `R-MILP -> MILP` warm-start approach (up to 10s for `R-MILP`, 300s total):
```bash
./ce_scheduler instances/automotive/inst1 300 -lb 10
```

## Input Instance Format

Lines starting with `#` or empty lines are ignored.

1. **Header line**:
   ```text
   <periodic_tasks> <sporadic_tasks> <processors> [minor_cycle]
   ```
   *(If `minor_cycle` is omitted, the frame size $F$ is computed automatically as the GCD of all periodic task periods, while the major cycle $P$ is their LCM).*

2. **Periodic tasks** (`<periodic_tasks>` lines, 0-indexed):
   ```text
   <period> <WCET>
   ```

3. **Sporadic tasks** (`<sporadic_tasks>` lines):
   ```text
   <min_inter_arrival> <WCET> <deadline> <jitter>
   ```

4. **Task relations** (remaining lines, between 0-indexed periodic tasks `task1` and `task2`):
   - `x <task1> <task2>`: Mutual exclusion relation ($\mathcal{E}$) between `task1` and `task2`.
   - `p <task1> <task2>`: Precedence relation ($\mathcal{P}$, `task1` precedes `task2`; both tasks must share the same period).

## Citation

If you use this code or benchmark instances in your research, please cite our paper:

```bibtex
@article{nogueira2026unified,
  title   = {Unified Scheduling of Periodic and Sporadic Tasks in Multicore Bare-Metal Systems via {MILP}},
  author  = {Nogueira, Bruno and Rocha, Pedro and Silva Junnior, Cleiber and Lima, Alfredo and Paes, Rodrigo and Silva, Francisco Airton},
  journal = {Journal of Systems Architecture},
  year    = {2026}
}
```
