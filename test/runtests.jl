using ParametrizedSurfaces
import ParametrizedSurfaces as PS
using Test
using StaticArrays
import Random: Xoshiro
using LinearAlgebra
using CoordinateTransformations
using Random
using Rotations

@testset "Sphere" begin
    radius = 10.0
    sphere = PS.Sphere(radius)
    rng = Xoshiro(2)
    θ = 2pi*rand(rng)
    ϕ = pi*rand(rng)
    uv = @SVector[θ,ϕ]
    pt = sphere(uv)
    @test surface_normal(sphere, uv) ≈ normalize(pt)
    @test principal_curvatures(sphere, uv) ≈ @SVector[1/radius, 1/radius]
    @test shape_operator(sphere, uv) ≈ PS.shape_operator_ad(sphere, uv)
end

@testset "paraboloid" begin
    function paraboloid(uv)
        u,v = uv
        z = 8*(1/2)*u^2 + 5*(1/2)*v^2 + 6*u*v -10
        @SVector[u,v,z]
    end
    @test shape_operator(paraboloid, [0,0]) ≈ -[8 6; 6 5]

    rng = Xoshiro(43)
    # uv = 0.1*randn(rng,2)
    uv = zeros(2)
    @test mean_curvature(paraboloid, uv) ≈ PS.mean_curvature_tr(paraboloid, uv)
    @test gaussian_curvature(paraboloid, uv) ≈ PS.gaussian_curvature_det(paraboloid, uv)

    @test shape_operator(paraboloid, uv) ≈ PS.shape_operator_ad(paraboloid, uv)

    K = gaussian_curvature(paraboloid, uv)
    H = mean_curvature(paraboloid, uv)
    k1, k2 = sort(principal_curvatures(paraboloid, uv))
    @test k1 ≈ H - √(H^2 -K)
    @test k2 ≈ H + √(H^2 -K)

    rot = rand(rng, AngleAxis)
    offset = SVector{3}(randn(rng, 3))
    trafo = AffineMap(rot, offset)
    @test principal_curvatures(paraboloid, uv) ≈ principal_curvatures(trafo∘paraboloid, uv)
    @test surface_normal(trafo∘paraboloid, uv) ≈ rot*surface_normal(paraboloid, uv)
end

@testset "Torus" begin
    @test PS.Torus(1,2)([0,0]) ≈ [3,0,0]
    @test PS.Torus(1,2)([0,pi]) ≈ [-3,0,0]
    @test PS.Torus(1,2)([pi,0]) ≈ [1,0,0]
    @test PS.Torus(1,2)([pi,pi]) ≈ [-1,0,0]
    r = 1.3
    R = 3.9
    t = PS.Torus(r,R)
    uv = randn(Xoshiro(345), 2)
    @test principal_curvatures(t, [0.0,0.0]) ≈ [-1/r, -1/(r+R),]
    @test principal_curvatures(t, [0.0,pi ]) ≈ [-1/r, -1/(r+R),]
    @test principal_curvatures(t, [pi ,0.0]) ≈ [-1/r, 1/(R-r),]   
    @test principal_curvatures(t, [pi ,pi])  ≈ [-1/r, 1/(R-r),]   
    @test shape_operator(t, uv) ≈ PS.shape_operator_ad(t, uv)
end

@testset "half sphere" begin
    function half_sphere10(uv)
        u,v = uv
        R=10.0
        z = sqrt(R^2 - u^2 - v^2)
        @SVector[u,v,z]
    end

    uv = @SVector[0.0,0]
    @inferred first_fundamental_form(half_sphere10, uv)
    @inferred second_fundamental_form(half_sphere10, uv)
    @inferred shape_operator(half_sphere10,uv)
    @inferred gaussian_curvature(half_sphere10, uv)
    @inferred mean_curvature(half_sphere10, uv)
    @inferred principal_curvatures(half_sphere10, uv)

    @test shape_operator(half_sphere10, uv) ≈ [0.1 0.0; 0.0 0.1]
    @test shape_operator(half_sphere10, uv) ≈ PS.shape_operator_ad(half_sphere10, uv)

    @test principal_curvatures(half_sphere10,uv) ≈ eigvals(shape_operator(half_sphere10, uv))

    rng = Xoshiro(23)
    @test shape_operator(half_sphere10, randn(rng, 2)) ≈ shape_operator(half_sphere10, zeros(2))
    # uv = rand(rng,2)

    uv = randn(rng, 2)
    @test eigvals(shape_operator(half_sphere10, uv)) ≈ principal_curvatures(half_sphere10, uv)
    @test shape_operator(half_sphere10, uv) ≈ PS.shape_operator_ad(half_sphere10, uv)
end

@testset "cone" begin
    function cone(uv)
        u,v = uv
        z = sqrt(u^2 + v^2)
        @SVector[u,v,z]
    end
    @test shape_operator(cone, [1,0]) ≈ [0 0; 0 -1/sqrt(2)]
    @test principal_curvatures(cone, [1,0]) ≈ [-1/sqrt(2), 0.0]
    rng = Xoshiro(42)
    for _ in 1:10
        uv = randn(rng, 2)
        r = norm(uv)
        @test eigvals(shape_operator(cone, uv)) ≈ [-inv(sqrt(2)*r), 0.0]
        m = shape_operator(cone, uv)
        @test m*uv ≈ [0, 0] atol=1e-15
        @test shape_operator(cone, uv) ≈ PS.shape_operator_ad(cone, uv)
    end
end
