import logging
import math
from ngsolve import *
import numpy as np

new_line = "\n"

def ref_element_vertices_val(gfu, vertices, gfu_index, element_type):
    val = []

    if element_type == "TRIG":
        # coefficients of the linear equation
        a = gfu.vec[gfu_index]
        b = gfu.vec[gfu_index + 1]
        c = gfu.vec[gfu_index + 2]
        for j in range(len(vertices)):
            val.append((a - b - c) + (3 * b + c) * vertices[j][0] + (2 * c) * vertices[j][1])

    if element_type == "TET":
        a = gfu.vec[gfu_index]
        b = gfu.vec[gfu_index + 1]
        c = gfu.vec[gfu_index + 2]
        d = gfu.vec[gfu_index + 3]
        for j in range(len(vertices)):
            eq = ((a - b - 2 * c - 4 * d) + (4 * b + 2 * c + 4 * d) * vertices[j][0] + (6 * c + 4 * d) * vertices[j][1] + 8 * d * vertices[j][2])
            val.append(eq)

    return val


def vertices_gfu_val(gfu, mesh, element_type="TRIG"):
    dof_per_element = int(len(gfu.vec) / mesh.ne)
    gfu_vertices_val = np.zeros((mesh.ne, dof_per_element))

    if element_type == "TRIG":
        vertices = [(0, 0), (1, 0), (0, 1)]
    if element_type == "TET":
        vertices = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]

    for i in range(mesh.ne):
        gfu_index = i * dof_per_element
        gfu_vertices_val[i, :] = ref_element_vertices_val(gfu, vertices, gfu_index, element_type)

    return gfu_vertices_val


def bound_limiter_Joshaghani(gfu, mesh, bounds):
    (r1, r2) = bounds
    averaged_val_gfu = gfu.vec.FV().NumPy()
    number_of_elements = mesh.ne
    dof_per_element = int(len(averaged_val_gfu) / number_of_elements)
    quad_val_gfu = vertices_gfu_val(gfu, mesh)

    quad_min_val = quad_val_gfu.min(axis=1)
    quad_max_val = quad_val_gfu.max(axis=1)

    theta = np.ones(number_of_elements, dtype=float)

    for i in range(number_of_elements):
        nn = dof_per_element * i
        if gfu.vec[nn] < r1:
            gfu.vec[nn] = r1
        if gfu.vec[nn] > r2:
            gfu.vec[nn] = r2
        if (quad_min_val[i] < r1):
            theta[i] = (gfu.vec[nn]-r1)/(gfu.vec[nn]-quad_min_val[i])
        if (quad_max_val[i] > r2):
            theta2 = (gfu.vec[nn]-r2)/(gfu.vec[nn]-quad_max_val[i])
            theta[i] = min(theta[i], theta2)
        if theta[i] < 1:
            for k in range(1, dof_per_element):
                gfu.vec[nn+k] = theta[i]*gfu.vec[nn+k]


def newton_solve(s, mesh, rhs, b, sold, a, As, mstar, w):
    rhs.data = b.mat * sold.vec
    rhs.data -= b.mat * s.vec
    a.Apply(s.vec, As)
    rhs.data -= As
    a.AssembleLinearization(s.vec)

    mstar.AsVector().data = b.mat.AsVector() + a.mat.AsVector()
    invmat = mstar.Inverse()
    w.data = invmat * rhs
    s.vec.data += w
    #bound_limiter_Joshaghani(s.components[0], mesh, (0, 1))
    wnorm = w.Norm()
    print("|w| = {:7.3e} ".format(wnorm), end="")
    return wnorm


def error(A, B, mesh):
    return sqrt(Integrate((A - B) ** 2, mesh))


def update_iterate(UN, UIter):
    UIter.vec.data = UN.vec.data


def timestep(UN, UOld):
    UOld.vec.data = UN.vec.data
