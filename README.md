# Numerical Methods for Bessel Functions

A MATLAB project exploring numerical approximation of Bessel functions of the first kind, with an emphasis on interpolation, quadrature, root-finding, and an application to circular membrane vibration.

## Features

- Polynomial interpolation using a Vandermonde system and the Lagrange form
- Interpolation using Chebyshev-distributed nodes
- Cubic spline interpolation, including a custom natural-spline implementation
- Numerical evaluation of Bessel functions using the integral representation
- Composite trapezoidal and Simpson's rules, plus Gauss–Legendre quadrature
- Root-finding with bisection, Newton–Raphson, the secant method, and MATLAB's `fzero`
- A membrane-vibration example that converts Bessel-function roots into natural frequencies

## Requirements

- MATLAB (tested conceptually with modern MATLAB releases; scripts use local functions only where supported)
- No third-party toolboxes are intended to be required

## Run

From MATLAB, set the current folder to the repository root and run:

```matlab
run('main.m')
```

The script creates comparison plots and prints root-finding results. To run only the membrane example, run:

```matlab
run('examples/membrane_vibration.m')
```

## Project structure

```text
bessel-numerical-methods/
├── main.m
├── README.md
├── LICENSE
├── .gitignore
├── src/
│   ├── besselIntegral.m
│   ├── besselTrapezoidal.m
│   ├── besselSimpson.m
│   ├── besselGaussLegendre.m
│   ├── legendreCoefficients.m
│   ├── interpolateVandermonde.m
│   ├── interpolateLagrange.m
│   ├── interpolateNaturalSpline.m
│   └── findBesselRoots.m
└── examples/
    └── membrane_vibration.m
```

## Numerical notes

- The tabulated values in `main.m` are rounded reference data. They are used to compare interpolation strategies, not as exact function values.
- High-degree polynomial interpolation on equally spaced nodes can be ill-conditioned. The Vandermonde method is included for comparison, not as a recommended general-purpose algorithm.
- The custom quadrature methods approximate the integral representation of the Bessel function. Accuracy depends on the number and distribution of quadrature nodes.
- Runtime measurements are machine- and load-dependent; compare them over repeated runs rather than treating a single timing as definitive.
- The membrane model assumes a uniform, circular membrane with constant tension and idealised boundary conditions. Real measurements can differ due to material properties, damping, and mounting conditions.

## Attribution

The tabulated reference values are based on *Handbook of Mathematical Functions* by Abramowitz and Stegun. The project implements numerical methods for learning and comparison; it is not intended as a validated engineering solver.
