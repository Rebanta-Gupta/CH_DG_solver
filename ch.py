from netgen.geom2d import unit_square
from ngsolve import *
import netgen.gui

from weakform import build_weakform
from solver import newton_solve
from amr import do_amr, save_vtk, REFINE_EVERY, MAX_NDOF
from ic_bc import set_initial_conditions
import user_settings as params


class CahnHilliardSolver:
    """CH/DG solver"""

    def __init__(self, mesh):
        self.mesh = mesh
        self.build_fes_and_forms()

    def build_fes_and_forms(self):
        """(Re-)build FE space and all dependent objects on the current mesh.
        Called once at construction and again after every mesh.Refine()."""
        V = L2(self.mesh, order=params.order, dgjumps=True)
        self.fes = FESpace([V, V])

        self.s    = GridFunction(self.fes)
        self.sold = GridFunction(self.fes)

        self.a, self.b = build_weakform(self.fes)
        self.mstar     = self.b.mat.CreateMatrix()

        self.space_flux = HDiv(self.mesh, order=params.order)
        self.gf_flux    = GridFunction(self.space_flux, "flux")

        self.rhs = self.s.vec.CreateVector()
        self.As  = self.s.vec.CreateVector()
        self.w   = self.s.vec.CreateVector()

    def _rebuild_solver_objects(self):
        """Rebuild a, b and work vectors without touching s / sold.
        Called after the IC is set so the forms see the populated solution."""
        self.a, self.b = build_weakform(self.fes)
        self.mstar     = self.b.mat.CreateMatrix()
        self.rhs       = self.s.vec.CreateVector()
        self.As        = self.s.vec.CreateVector()
        self.w         = self.s.vec.CreateVector()

    def time_step(self):
        """One implicit Euler step via Newton iteration."""
        self.sold.vec.data = self.s.vec.data
        wnorm = 1e99
        if params.solver == "Newton":
            while wnorm > 20:
                wnorm = newton_solve(self.s, self.mesh, self.rhs, self.b,
                                     self.sold, self.a, self.As, self.mstar,
                                     self.w)


def main():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.04))
    ch   = CahnHilliardSolver(mesh)

    set_initial_conditions(ch)
    ch._rebuild_solver_objects()

    Draw(ch.s.components[0], ch.mesh, "c")
    Draw(ch.s.components[1], ch.mesh, "mu")

    t    = 0.0
    step = 0
    save_vtk(ch.mesh, ch.s, step)

    while t < params.tend:
        print("\n\nt = {:10.6e}  (step {:d})".format(t, step))

        ch.time_step()

        t    += params.tau
        step += 1
        Redraw(blocking=False)
        save_vtk(ch.mesh, ch.s, step)

        if (params.amr
                and REFINE_EVERY > 0
                and step % REFINE_EVERY == 0
                and ch.fes.ndof < MAX_NDOF):
            do_amr(ch)
            Draw(ch.s.components[0], ch.mesh, "c")
            Draw(ch.s.components[1], ch.mesh, "mu")


if __name__ == "__main__":
    main()
