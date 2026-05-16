import numpy as np
from ngsolve import CoefficientFunction


def set_initial_conditions(solver):
    """ Each element independently gets 0 or 1.

    Parameters
    ----------
    solver : CahnHilliardSolver  (from ch.py)
    """
    gf_c = solver.s.components[0]
    vec  = gf_c.vec.FV().NumPy()

    n_el       = solver.mesh.ne
    dof_per_el = len(vec) // n_el

    rng    = np.random.default_rng(seed=42)
    values = rng.integers(0, 2, size=n_el).astype(float)

    for i in range(n_el):
        vec[i * dof_per_el] = values[i]
        for k in range(1, dof_per_el):
            vec[i * dof_per_el + k] = 0.0

    solver.s.components[1].Set(CoefficientFunction(0.0))
    solver.sold.vec.data = solver.s.vec.data

    print("IC set: random per-element binary ({:d} elements)".format(n_el))
