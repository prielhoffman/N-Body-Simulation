# MPI N-Body Simulation

A **parallel N-body simulation in C using MPI**, modeling 2D gravitational interactions between hundreds of bodies and measuring the performance of distributed computation across multiple processes.

This project demonstrates **parallel programming**, **numerical simulation**, and **performance analysis** in a high-computation scientific workload.

## Tech Stack

- **Language:** C
- **Parallel Framework:** MPI
- **Concepts:** parallel computing, distributed computation, simulation, numerical methods, performance analysis

## What the Project Does

The project simulates the gravitational interaction of **992 bodies** in a 2D domain.

Each body affects every other body according to Newtonian gravity, and the system evolves over time by repeatedly:

- computing gravitational forces
- updating velocities
- updating positions
- handling bodies that leave the simulation domain

The force calculation follows the naive **O(n²)** approach, making it computationally expensive and a strong candidate for parallelization with MPI. :contentReference[oaicite:1]{index=1}

## Main Features

- parallel N-body simulation in C
- MPI-based distribution of computation
- 2D gravitational interaction model
- random initialization of positions and velocities
- iterative time-step simulation
- handling of bodies leaving the simulation domain
- performance measurement for different process counts
- analysis of speedup and efficiency

## Simulation Model

The simulation models a set of bodies moving under mutual gravitational attraction.

### Initial Conditions

The project initializes:

- **992 stars**
- positions inside a 2D square domain
- random velocities within a defined range

The system then evolves over a sequence of time steps using double-precision calculations. :contentReference[oaicite:2]{index=2}

### Force Computation

For each body, the program computes the total gravitational effect of all other bodies.

Because every body interacts with every other body, the baseline computation cost is **quadratic** in the number of bodies.

This makes the simulation a useful benchmark for parallel execution. :contentReference[oaicite:3]{index=3}

### State Update

After computing forces, the simulation updates:

- body velocities
- body positions

This process is repeated for multiple time steps to simulate the system over time. :contentReference[oaicite:4]{index=4}

### Domain Handling

The implementation includes logic for bodies that move outside the simulation domain, such as reflection or wrapping behavior, depending on the chosen strategy. :contentReference[oaicite:5]{index=5}

## Parallelization

The project uses **MPI** to distribute the computation across multiple processes.

Parallel execution allows the expensive force calculations to be shared between processes, making it possible to reduce runtime and study the scalability of the simulation.

This project is especially useful for understanding how compute-heavy numerical problems behave under distributed execution.

## What I Implemented

This project focused on combining scientific simulation with parallel performance analysis.

Key implementation areas included:

- modeling 2D gravitational interactions
- implementing iterative state updates
- parallelizing the computation with MPI
- handling initialization of simulation state
- using double precision for numerical calculations
- measuring runtime under different process counts
- analyzing speedup and efficiency of the parallel program

## Why This Project Matters

This project demonstrates practical understanding of:

- MPI-based parallel programming
- distributed computation for compute-intensive workloads
- numerical simulation in C
- the performance impact of O(n²) algorithms
- measuring scalability, speedup, and efficiency
- structuring scientific computation for parallel execution

It is a strong technical project because it combines low-level C programming, parallel systems concepts, and quantitative performance evaluation.

## Performance Analysis

The project measures execution time for different numbers of MPI processes and uses the results to evaluate:

- runtime improvement
- speedup
- efficiency

This turns the project from a pure simulation into a performance-oriented parallel computing exercise. :contentReference[oaicite:6]{index=6}

## How to Run

Compile the program:

```bash
mpicc -o nbody nbody.c
```

Run the simulation:

```bash
mpirun -np <number_of_processes> ./nbody
```

Example:

```bash
mpirun -np 4 ./nbody
```

Make sure MPI is installed and configured correctly in your environment before compiling and running the program. :contentReference[oaicite:7]{index=7}

## Core Concepts Practiced

- MPI
- parallel computing
- distributed processing
- scientific simulation
- N-body modeling
- performance measurement
- speedup and efficiency analysis
- C programming for high-compute workloads

## Key Takeaways

Through this project, I strengthened my understanding of:

- how to parallelize compute-heavy simulations
- how O(n²) workloads behave under distributed execution
- how to evaluate scalability using runtime, speedup, and efficiency
- how to structure numerical simulation code in C
- how MPI can be used to coordinate parallel scientific computation

## Future Improvements

Possible next steps for the project:

- compare naive and optimized force-computation strategies
- visualize body movement over time
- log performance results more systematically
- test larger body counts and longer simulations
- compare scaling behavior across more process counts
- add charts summarizing speedup and efficiency results
