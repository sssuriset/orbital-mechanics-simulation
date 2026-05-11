# Energy Conservation and Timestep Sensitivity in Numerical Orbital Simulations

This project compares three numerical integration methods for simulating two-body orbital motion around Earth: Euler, semi-implicit Euler, and Velocity Verlet. The goal is to evaluate how timestep size affects orbital stability, energy conservation, angular momentum conservation, and orbital period accuracy.

## Methods

The simulation models a satellite in Earth orbit using Newtonian gravity. Each method is tested with timestep values of 1, 5, 10, 30, 60, 120, and 300 seconds over a seven-day simulation. The same initial position and velocity are used for all methods.

The tested methods are:
- Euler
- Semi-implicit Euler
- Velocity Verlet

## Metrics

The simulation records:
- total mechanical energy drift
- angular momentum drift
- final orbital radius error
- orbital period error
- runtime

## Results

The results show that timestep size strongly affects numerical stability. Euler integration accumulates the largest energy drift as timestep increases. Semi-implicit Euler improves stability but still shows significant error at larger timesteps. Velocity Verlet maintains much lower energy drift across the tested timestep range and produces more stable orbital trajectories.

## Output Files

The script saves result figures and a CSV file in the `results` folder.

Generated files include:
- `orbitcomparison.png`
- `energycomparison.png`
- `angmomcomparison.png`
- `timestepsensitivity.png`
- `perioderror.png`
- `results.csv`
