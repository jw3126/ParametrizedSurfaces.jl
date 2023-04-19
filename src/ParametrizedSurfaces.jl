module ParametrizedSurfaces
using StaticArrays
using ForwardDiff
using ArgCheck
using LinearAlgebra

export first_fundamental_form, second_fundamental_form
export shape_operator
export principal_curvatures, gaussian_curvature, mean_curvature
export surface_normal

function calc_EFG(f,uv)
    uv = float(uv)
    J = ForwardDiff.jacobian(f, uv)
    f_u = J[:,1]
    f_v = J[:,2]
    E = f_u⋅f_u
    F = f_u⋅f_v
    G = f_v⋅f_v
    (;E,F,G)
end

function first_fundamental_form(f, uv)
    (;E,F,G) = calc_EFG(f, uv)
    return @SMatrix[E F; F G]
end

function calc_LMN(f, uv)
    m = second_fundamental_form(f,uv)
    L = m[1,1]
    M = m[2,1]
    N = m[2,2]
    return (;L,M,N)
end

function second_fundamental_form(f, uv)
    uv = float(uv)
    n = surface_normal(f, uv)
    h = ForwardDiff.hessian(uv->f(uv)⋅n, uv)
    return h
end

function shape_operator(f, uv)
    F1 = first_fundamental_form(f, uv)
    F2 = second_fundamental_form(f, uv)
    return -F2*inv(F1)
end

function shape_operator_ad(f,uv)
    J = ForwardDiff.jacobian(f, uv)
    Jn = ForwardDiff.jacobian(uv) do uv_
        surface_normal(f, uv_)
    end
    J \ Jn
end

function gaussian_curvature(f,uv)
    (;L,M,N) = calc_LMN(f,uv)
    (;E,F,G) = calc_EFG(f,uv)
    K = (L*N - M^2) / (E*G - F^2)
    return K
end

function gaussian_curvature_det(f, uv)
    det(shape_operator(f, uv))
end

function mean_curvature(f,uv)
    (;L,M,N) = calc_LMN(f,uv)
    (;E,F,G) = calc_EFG(f,uv)
    H = (G*L-2*F*M + E*N) / (2*(E*G-F^2))
    return -H
end

function mean_curvature_tr(f, uv)
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

function cross_uv(f, uv)
    uv = float(uv)
    J = ForwardDiff.jacobian(f, uv)
    f_u = J[:,1]
    f_v = J[:,2]
    return cross(f_u, f_v)
end

function surface_normal(f, uv)
    n = cross_uv(f, uv)
    return normalize(n)
end

#################################################################################
#### Example surfaces
################################################################################
struct Torus
    r::Float64
    R::Float64
end
function (torus::Torus)(uv)
    (;r,R) = torus
    γ,Γ = uv
    sinγ,cosγ = sincos(γ)
    sinΓ,cosΓ = sincos(Γ)
    x=(R+r*cosγ)*cosΓ
    y=(R+r*cosγ)*sinΓ
    z=r*sinγ
    @SVector[x,y,z]
end

struct Sphere
    radius::Float64
end
function (sphere::Sphere)(uv)
    (;radius) = sphere
    θ,ϕ = uv
    x = radius*sin(θ)*cos(ϕ)
    y = radius*sin(θ)*sin(ϕ)
    z = radius*cos(θ)
    @SVector[x,y,z]
end

# TODO
# struct Cone
# end

# struct Cylinder
# end

end#module
