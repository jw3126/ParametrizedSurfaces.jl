using ParametrizedSurfaces
using Test
using StaticArrays
import Random: Xoshiro
using LinearAlgebra

@testset "sphere" begin
    function sphere10(uv)
        u,v = uv
        R=10.0
        z = sqrt(R^2 - u^2 - v^2)
        @SVector[u,v,z]
    end

    uv = @SVector[0.0,0]
    @inferred first_fundamental_form(sphere10, uv)
    @inferred second_fundamental_form(sphere10, uv)
    @inferred shape_operator(sphere10,uv)
    @inferred gaussian_curvature(sphere10, uv)
    @inferred mean_curvature(sphere10, uv)
    @inferred principal_curvatures(sphere10, uv)

    @test shape_operator(sphere10, uv) ≈ [-0.1 0.0; 0.0 -0.1] # TODO sign correct?

    @test principal_curvatures(sphere10,uv) ≈ eigvals(shape_operator(sphere10, uv))

    rng = Xoshiro(23)
    @test shape_operator(sphere10, randn(rng, 2)) ≈ shape_operator(sphere10, zeros(2))

    uv = randn(rng, 2)
    @test eigvals(shape_operator(sphere10, uv)) ≈ principal_curvatures(sphere10, uv)
end

@testset "cone" begin
    function cone(uv)
        u,v = uv
        z = sqrt(u^2 + v^2)
        @SVector[u,v,z]
    end

    @test shape_operator(cone, [1,0]) ≈ [0 0; 0 1/sqrt(2)]
    @test principal_curvatures(cone, [1,0]) ≈ [0, 1/sqrt(2)]

    rng = Xoshiro(42)
    for _ in 1:10
        uv = randn(rng, 2)
        r = norm(uv)
        @test eigvals(shape_operator(cone, uv)) ≈ [0, inv(sqrt(2)*r)]
        m = shape_operator(cone, uv)
        @test m*uv ≈ [0, 0] atol=1e-15
    end
end

@testset "paraboloid" begin
    function paraboloid(uv)
        u,v = uv
        z = 8*(1/2)*u^2 + 5*(1/2)*v^2 + 6*u*v -10
        @SVector[u,v,z]
    end
    @test shape_operator(paraboloid, [0,0]) ≈ [8 6; 6 5]
end

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

@testset "torus" begin
    @test Torus(1,2)([0,0]) ≈ [3,0,0]
    @test Torus(1,2)([0,pi]) ≈ [-3,0,0]
    @test Torus(1,2)([pi,0]) ≈ [1,0,0]
    @test Torus(1,2)([pi,pi]) ≈ [-1,0,0]
    r = 1.2
    R = 2.8
    t = Torus(r,R)
    @test principal_curvatures(t, [0.0,0.0]) ≈ [1/(r+R), 1/r]
    @test principal_curvatures(t, [0.0,pi ]) ≈ [1/(r+R), 1/r]
    @test principal_curvatures(t, [pi ,0.0]) ≈ [-1/(R-r), 1/r]   
    @test principal_curvatures(t, [pi ,pi])  ≈ [-1/(R-r), 1/r]   
end
