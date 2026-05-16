# Cahn–Hilliard DG Solver

Solves the Cahn–Hilliard equation using a discontinuous Galerkin (DG) method with adaptive mesh refinement, implemented in [NGSolve](https://ngsolve.org).

---

## Equations

The Cahn–Hilliard equation describes phase separation of a binary mixture. The phase field $c \in [0,1]$ and chemical potential $\mu$ satisfy:

$$\frac{\partial c}{\partial t} + \nabla \cdot (-M \nabla \mu) = 0$$

$$\mu = \frac{dF}{dc} - \lambda \nabla^2 c$$

where the bulk free energy density is the double-well potential:

$$F(c) = 100\, c^2 (1 - c)^2$$

giving:

$$\frac{dF}{dc} = 200\, c\,(1 - c)\,(1 - 2c)$$

---

## Method

- **Spatial discretisation:** Interior Penalty DG
- **Time integration:** Implicit Euler with Newton iteration
- **AMR:** Zienkiewicz–Zhu gradient-recovery error estimator on $\nabla c$

---

## Sample Output

![Phase field evolution](anim.gif)
