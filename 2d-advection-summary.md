# 2D Advection Simulation Program

## Program Description

This C program simulates 2D advection of a Gaussian distribution u(x,y) at a fixed velocity, with modifications to include atmospheric boundary layer (ABL) conditions.

## Compilation and Execution

Compile: `gcc -o advection2D -std=c99 advection2D.c -lm`

## Key Components

### Constants and Parameters

- Grid: 1000x1000 points
- Domain: x [0, 30], y [0, 30]
- Gaussian initial conditions: centered at (3, 15) with widths (1, 5)
- Time stepping: CFL = 0.9, 1000 steps
- ABL parameters: friction velocity = 0.1, von Karman constant = 0.41, roughness length = 1.0

### Main Functions and Algorithms

1. Grid initialization
2. Setting up Gaussian initial conditions
3. Time-stepping loop:
   - Apply boundary conditions
   - Calculate rate of change (dudt) using leftward difference
   - Update u values
4. Output results to files

### OpenMP Parallelization

Several loops are parallelized using OpenMP:
- Grid initialization
- Initial condition setup
- Boundary condition application
- Rate of change calculation

### File Output

- `initial.dat`: Initial u(x,y) values
- `final.dat`: Final u(x,y) values
- `avg_file.dat`: Average u values along x-axis

## Key Features

1. Uses CFL condition for time step calculation
2. Implements atmospheric boundary layer (ABL) conditions
3. Calculates and outputs average u values
4. Utilizes OpenMP for parallel computation

## Notes and Potential Improvements

1. The velocity in x-direction (velx) is updated based on ABL conditions
2. Some loops are not parallelized to maintain specific output order
3. Consider using more efficient I/O methods for large datasets
4. Potential for further optimization and parallelization

