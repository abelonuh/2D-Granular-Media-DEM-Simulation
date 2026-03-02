# 2D Granular Media Simulation Using the Discrete Element Method (DEM)

## Overview

This project implements a two-dimensional Discrete Element Method (DEM) solver in Fortran 90 to simulate granular particle interactions under gravity.

The work was completed collaboratively by **Antoine Millet–Cot** and **Abel Inalegwu Onuh** as part of a numerical simulation course at Mines Saint-Étienne (January 2026).

The solver models particle–particle and particle–wall collisions using a linear spring–dashpot contact model and integrates particle motion using an explicit time integration scheme.

---

## Physical Model

Each particle obeys Newton’s second law:

F_net = F_gravity + F_normal + F_damping

Contact mechanics are modeled using a soft-sphere linear spring–dashpot formulation:

- Normal elastic force: F = k_n δ
- Viscous damping force: F = c_n δ_dot

The analytical damped harmonic oscillator solution was used to validate the collision model.

---

## Numerical Implementation

- Language: Fortran 90
- Collision detection optimized to avoid redundant pair calculations
- Explicit Euler time integration
- Stability condition enforced: Δt << t_c

The theoretical contact duration was computed and used to justify the selected timestep.

---

## Validation

Validation was performed using a head-on collision benchmark:

- Numerical overlap evolution matched analytical solution
- Stability analysis confirmed divergence for large Δt
- Energy dissipation verified in damped simulations

---

## Results

The solver successfully simulates:

- Particle settling under gravity
- Inelastic collisions with energy dissipation
- Multi-body granular interactions

See the `results/` folder for animation output.

---

## Compilation

Use:

gfortran collisions.f90 vfinw.f90 -o main

Run:

./main

---

## Authors

Antoine Millet–Cot  
Abel Inalegwu Onuh  

Mines Saint-Étienne, 2026
