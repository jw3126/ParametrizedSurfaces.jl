module ParametrizedSurfaces
using StaticArrays
using ForwardDiff
using ArgCheck
using LinearAlgebra

export first_fundamental_form, second_fundamental_form
export shape_operator
export principal_curvatures, gaussian_curvature, mean_curvature

function first_fundamental_form(f, uv)
    uv = float(uv)
    J = ForwardDiff.jacobian(f, uv)
    f_u = J[:,1]
    f_v = J[:,2]
    E = f_u⋅f_u
    F = f_u⋅f_v
    G = f_v⋅f_v
    @SMatrix[E F; F G]
end

function second_fundamental_form(f, uv)
    uv = float(uv)
    J = ForwardDiff.jacobian(f, uv)
    f_u = J[:,1]
    f_v = J[:,2]
    n = normalize(cross(f_u, f_v))
    h = ForwardDiff.hessian(uv->f(uv)⋅n, uv)
    return h
end

function shape_operator(f, uv)
    F1 = first_fundamental_form(f, uv)
    F2 = second_fundamental_form(f, uv)
    F2*inv(F1)
end

function gaussian_curvature(f, uv)
    det(shape_operator(f, uv))
end

function mean_curvature(f, uv)
    tr(shape_operator(f,uv))/2
end

function eigenvalues2x2(M)
    m = tr(M)/2
    d = det(M)
    Δ = sqrt(m^2 - d)
    @SVector[m - Δ, m + Δ]
end

function principal_curvatures(f, uv)
    M = shape_operator(f,uv)
    eigenvalues2x2(M)
end

end#module
