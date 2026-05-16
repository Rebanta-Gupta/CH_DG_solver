"""Adaptive mesh refinement utilities for the CH DG solver.

Provides:
  - calc_error_and_mark : ZZ gradient-recovery error estimator
  - do_amr              : refine mesh and prolong solution
  - save_vtk            : write one VTU snapshot per step
"""

from ngsolve import *
import user_settings as params

# ── AMR tuning ────────────────────────────────────────────────────────────────
REFINE_EVERY    = 5       # run AMR every this many time steps (0 = disabled)
REFINE_FRACTION = 0.25    # mark elements with local error > fraction * max_error
MAX_NDOF        = 80_000  # stop refining beyond this DOF count


def calc_error_and_mark(mesh, fes, gf_c, space_flux, gf_flux):
    """Zienkiewicz-Zhu gradient-recovery error estimator on the phase field c.

    Recovers grad(c) into an H(div) space, computes the element-wise squared
    difference, and marks elements whose local error exceeds
    REFINE_FRACTION * max_error for refinement.

    Parameters
    ----------
    mesh       : ngsolve.Mesh
    fes        : ngsolve.FESpace
    gf_c       : GridFunction  (phase field component c)
    space_flux : HDiv space
    gf_flux    : GridFunction on space_flux

    Returns
    -------
    global_err : float
    """
    space_flux.Update()
    gf_flux.Update()

    gf_flux.Set(grad(gf_c))

    err   = (grad(gf_c) - gf_flux) * (grad(gf_c) - gf_flux)
    elerr = Integrate(err, mesh, VOL, element_wise=True)

    maxerr    = max(elerr)
    global_err = sqrt(sum(elerr))

    print("  [AMR] ndof = {:6d}  |  est. error = {:.3e}  |  max elem = {:.3e}"
          .format(fes.ndof, global_err, maxerr))

    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, elerr[el.nr] > REFINE_FRACTION * maxerr)

    return global_err


def do_amr(solver):
    """Run the ZZ estimator, refine marked elements, and prolong the solution.

    Parameters
    ----------
    solver : CahnHilliardSolver  (from ch.py)
    """
    calc_error_and_mark(
        solver.mesh, solver.fes,
        solver.s.components[0],
        solver.space_flux, solver.gf_flux,
    )

    # Keep strong references to old GridFunctions across the rebuild so that
    # Set(...) can interpolate from them after build_fes_and_forms() replaces
    # self.s and self.sold with new (empty) objects.
    s_old    = solver.s
    sold_old = solver.sold

    solver.mesh.Refine()
    solver.build_fes_and_forms()

    solver.s.components[0].Set(s_old.components[0])
    solver.s.components[1].Set(s_old.components[1])
    solver.sold.vec.data = solver.s.vec.data
    # s_old / sold_old go out of scope here, after Set() has finished.

    print("  [AMR] mesh refined  ->  new ndof = {:d}".format(solver.fes.ndof))


def save_vtk(mesh, s, step, subdivision=1):
    """Write one VTU snapshot with the step number in the filename.

    A fresh VTKOutput is created each call so NGSolve's internal counter
    is always 0 — filenames never collide even after mesh refinement.

    Parameters
    ----------
    mesh        : ngsolve.Mesh
    s           : GridFunction  (mixed c / mu)
    step        : int
    subdivision : int
    """
    if params.vtkoutput:
        fname = "solution/cahnhilliard_step{:05d}".format(step)
        VTKOutput(ma=mesh,
                  coefs=[s.components[1], s.components[0]],
                  names=["mu", "c"],
                  filename=fname,
                  subdivision=subdivision).Do()
