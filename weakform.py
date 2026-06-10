from ngsolve import *
from helpers import Jump, gradavg
from config import params


def build_weakform(fes):
    """Build the nonlinear CH bilinear form (a) and mass form (b).

    Parameters
    ----------
    fes : ngsolve.FESpace  (mixed L2 x L2 DG)

    Returns
    -------
    a : BilinearForm   nonlinear residual / Jacobian (not assembled)
    b : BilinearForm   mass-like term (assembled)
    """
    c, mu = fes.TrialFunction()
    q, v  = fes.TestFunction()
    n     = specialcf.normal(2)
    h     = specialcf.mesh_size
    alpha = 4 * params.order**2 / h

    tau   = params.tau
    M     = params.M
    lam   = params.lamdba

    # ── chemical potential flux ───────────────────────────────────────────────
    a = BilinearForm(fes)
    a += tau * M * grad(mu) * grad(q) * dx
    a += -tau * M * n * gradavg(mu) * Jump(q) * dx(skeleton=True)
    a += -tau * M * n * gradavg(q)  * Jump(mu) * dx(skeleton=True)
    a += tau * M * alpha * Jump(q) * Jump(mu) * dx(skeleton=True)

    # ── mu = dF/dc  (double-well F = 100 c^2 (1-c)^2) ────────────────────────
    a += mu * v * dx
    a += -200 * (c - 3 * c**2 + 2 * c**3) * v * dx

    # ── interface energy (Laplacian of c, DG) ────────────────────────────────
    a += -lam * grad(c) * grad(v) * dx
    a += lam * n * gradavg(c) * Jump(v) * dx(skeleton=True)
    a += lam * n * gradavg(v) * Jump(c) * dx(skeleton=True)
    a += -lam * alpha * Jump(c) * Jump(v) * dx(skeleton=True)

    # ── convection ────────────────────────────────────────────────────────────
    vel = CoefficientFunction((2, 2))
    a += tau * -c * vel * grad(q) * dx
    a += tau * c * IfPos(vel * n, vel * n, 0) * Jump(q) * dx(skeleton=True)
    a += tau * c * IfPos(vel * n, 0, vel * n) * Jump(q) * dx(skeleton=True)

    # ── mass matrix ───────────────────────────────────────────────────────────
    b = BilinearForm(fes)
    b += SymbolicBFI(c * q)
    b.Assemble()

    return a, b
