using QuasiMonteCarlo
using DelaunayTriangulation
using CairoMakie
using StableRNGs

begin
    idx = rand(1:10000)
    rng = StableRNG(idx)
    lb = (-π, -π)
    ub = (1*π, 1*π)
    points = QuasiMonteCarlo.sample(500, lb, ub, HaltonSample())
    ch = convex_hull(points)
        
    tri = triangulate(points; boundary_nodes=ch.vertices, rng)

    ch_points = [get_point(tri, i) for i in DelaunayTriangulation.get_vertices(ch)]
    fig1, ax, sc = lines(ch_points, color = :red, linewidth = 4)
    scatter!(ax, points)
    
    fig2 = triplot(tri)
    fig1
end

function clean_points(points; tol=1e-10)
    kept = Vector{Int}()
    for j in axes(points,2)
        x = points[:,j]
        if all(norm(x - points[:,i]) > tol for i in kept)
            push!(kept,j)
        end
    end
    return points[:,kept]
end

points = clean_points(points; tol=1e-8)


using DelaunayTriangulation

tri = triangulate(points)
hull = convex_hull(tri)

using Gmsh

gmsh.initialize()
gmsh.model.add("cloud_hull")

lc = 0.45   # target mesh size

pids = Int[]

for i in hull.vertices[1:end-1]
    x = points[1, i]
    y = points[2, i]
    push!(pids, gmsh.model.geo.addPoint(x, y, 0.0, lc))
end

lids = Int[]
n = length(pids)

for i in 1:n
    p1 = pids[i]
    p2 = pids[mod1(i+1, n)]
    push!(lids, gmsh.model.geo.addLine(p1, p2))
end

loop = gmsh.model.geo.addCurveLoop(lids)
surf = gmsh.model.geo.addPlaneSurface([loop])

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [surf], 1)
gmsh.model.setPhysicalName(2, 1, "domain")

gmsh.model.addPhysicalGroup(1, lids, 2)
gmsh.model.setPhysicalName(1, 2, "boundary")

gmsh.model.mesh.generate(2)
gmsh.write("cloud_hull.msh")
gmsh.finalize()

gmsh.initialize()
gmsh.open("cloud_hull.msh")
gmsh.fltk.run()   # launches GUI
gmsh.finalize()

using Gridap
using GridapGmsh

model = GmshDiscreteModel("cloud_hull.msh")

order = 3
reffe = ReferenceFE(lagrangian, Float64, order)

V = TestFESpace(
    model,
    reffe;
    dirichlet_tags = ["boundary"]
)

g(x) = 1
U = TrialFESpace(V, g)

degree = 2 * order
Ω = GridapGmsh.Triangulation(model)
dΩ = Measure(Ω, degree)

# Right-hand side
f(x) = 1.0

a(u, v) = ∫( ∇(v) ⋅ ∇(u) )dΩ
l(v)    = ∫( v * f )dΩ

op = AffineFEOperator(a, l, U, V)
A = Matrix(get_matrix(op))
b = get_vector(op)
u = A \ b
uh = solve(op)


using CairoMakie
using GridapMakie
using GLMakie

fig, ax, plt = plot(Ω, uh, shading=true)
Colorbar(fig[1,2], plt)