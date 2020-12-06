# SUSY-algebra
Python software to perform algebraic computation in supersymmetric quantum field theory in 4d

The main library is SUSYalgebra.py. It computes the Lagrangian in components from its definition in terms of superfields and covariant derivatives, as well as the equations of motion, the supercurrent, and the energy-momentum tensor. The output is either human-readable, or in the form of code that can be interpreted by Mathematica.

The following files are executable:
- Matter.py: free chiral multiplet
- Matter_local.py: free chiral multiplet with a multiplicative coupling superfield
- Gauge.py: pure gauge theory
- Gauge_local.py: pure gauge theory with a multiplicative coupling superfield
- SYM.py: N=2 super-Yang-Mills Lagrangian, including gauge and matter superfields
- SYM_local.py: N=2 super-Yang-Mills Lagrangian with a multiplicative coupling superfield
- Higher_deriv_1.py and Higher_deriv_2.py: Lagrangians consisting in derivatives of a dimensionless superfield
- SYM_complete.py: N=2 super-Yang-Mills theory with a multiplicative coupling superfield
- Tests.py: various tests of the package SUSYalgebra.py

This code has been used to compute the supercurrent and the energy-momentum tensor in N=2 supersymmetric Yang-Mills theory with a space-dependent cutoff. The results have been published at https://arxiv.org/abs/1707.05325

---

Note 1: This package needs the SymPy library for symbolic calculations (www.sympy.org)

Note 2: The superfields are represented by vectors, and the derivative operators in superspace act as matrix multiplication. This is a `brute force' method that could certainly be optimized by a more clever use of the supersymmetric algebra, but it is sufficient for many purposes.
