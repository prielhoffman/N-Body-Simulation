# N-Body Simulation Using MPI

This repository contains a parallelized C program for simulating a 2D N-body problem using MPI (Message Passing Interface). The simulation models the gravitational interactions of 992 stars in a square domain, with the computation of forces and movements performed in parallel across multiple processes.

## Overview

The program simulates the gravitational interactions of stars in a 2D space, using the naive O(n²) approach to compute forces between all pairs of stars. Initial positions and velocities are generated randomly, and the simulation runs for a specified number of time steps. Performance measurements are conducted to analyze the efficiency and effectiveness of parallel execution using MPI.

## Key Features

- **Parallel Computation**: Utilizes MPI to parallelize the simulation of gravitational interactions among stars.
- **Initial Conditions**: Randomly generates positions and velocities for 992 stars within a defined square domain.
- **Gravitational Calculation**: Computes gravitational forces and updates star positions based on Newtonian physics.
- **Domain Handling**: Includes a strategy for handling stars exiting the simulation domain (e.g., reflecting back or wrapping around).
- **Performance Measurement**: Assesses the performance of the parallel program on “hobbit” nodes with varying numbers of processes.

## How It Works

1. **Initial Setup**:
    - **Positions**: Randomly generate initial positions for 992 stars within a 100 light-year by 100 light-year domain.
    - **Velocities**: Assign initial velocities between 0.5v and 1.5v, where v is the average speed (200 km/sec), with uniform direction distributions.

2. **Simulation**:
    - **Gravitational Force Calculation**: For each star, compute gravitational forces exerted by every other star using Newton’s law of gravitation.
    - **Position Update**: Update positions and velocities of stars based on the computed forces.
    - **Domain Handling**: Implement a strategy for stars exiting the domain (e.g., reflection or wrapping).

3. **Execution**:
    - **Time Steps**: Choose an appropriate time step and total number of time steps to ensure the simulation runs for more than 2 minutes.
    - **Double Precision**: Use double precision for physical variables to ensure accurate computations.

4. **Performance Measurement**:
    - **Execution Time**: Measure the execution time using MPI functions on the VM platform.
    - **Analysis**: Evaluate performance metrics such as speedup and efficiency for different numbers of MPI processes.

## Running the Program

### 1. Compile the Program:

```bash
mpicc -o nbody nbody.c
```

### 2. Execute the Program with different numbers of MPI processes:

```bash
mpirun -np [number_of_processes] ./nbody
```
Example for 4 processes:
```bash
mpirun -np 4 ./nbody
```

## Notes
* Ensure that MPI is properly installed and configured on your system.
* The choice of time step and total number of steps should ensure a simulation duration of over 2 minutes.
* Decide on an appropriate strategy for handling stars that exit the domain.
