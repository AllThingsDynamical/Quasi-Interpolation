using QuasiMonteCarlo
using DelaunayTriangulation
using CairoMakie
using StableRNGs

begin
    idx = rand(1:10000)
    rng = StableRNG(idx)
    lb = (-π, -π)
    ub = (1*π, 1*π)
    points = QuasiMonteCarlo.sample(1_000, lb, ub, LatinHypercubeSample())
    tri = triangulate(points; rng)

    ch = convex_hull(points)
    ch_points = [get_point(tri, i) for i in DelaunayTriangulation.get_vertices(ch)]
    fig1, ax, sc = lines(ch_points, color = :red, linewidth = 4)
    scatter!(ax, points)
    
    fig2 = triplot(tri)
    fig1
end

function target_function(X::Matrix)
    M = size(X,2)
    Y = zeros(1,M)
    f = (x,y) -> sin(x)*sin(y) + sin(4*x)*sin(4*y)
    for i=1:M
        x = X[:,i]
        Y[1,i] = f(x[1], x[2])
    end
    return Y
end

function get_hull_points(points)
    ch = convex_hull(points)
    ch_points = [collect(get_point(tri, i)) for i in DelaunayTriangulation.get_vertices(ch)]
    return reduce(hcat, ch_points)
end

X = points
Y = target_function(X)

X_hull = get_hull_points(points)
Y_hull = target_function(X_hull)

