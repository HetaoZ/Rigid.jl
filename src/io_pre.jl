"""
Generate a rigid structure.
"""
function create_rigid(dim::Int, rho::Real, I::Real, xs::Vector...)
    s = RigidStructure(dim, length(xs[1]))
    
    for i in eachindex(s.nodes)
        s.nodes[i] = Node(i, [xs[k][i] for k in 1:dim], zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim))
    end
    
    s.system.x0 = [mean(xs[k]) for k in 1:dim]
    s.system.x = [mean(xs[k]) for k in 1:dim]
    s.system.theta = 0.
    s.system.omega = 0.
    s.system.domega = 0.

    if dim == 1
        s.boundary = [Convex(1,[1],[1.0]), Convex(2,[2],[1.0])]
        if s.nodes[1].x0[1] < s.nodes[2].x0[1]
            s.boundary[1].normal *= -1
        else
            s.boundary[2].normal *= -1
        end
    elseif dim == 2
        for i in eachindex(s.boundary)
            s.boundary[i] = Convex(i, [i,i+1], [0.,0.])
        end
        s.boundary[end].link[2] = 1
        for i in eachindex(s.boundary)
            s.boundary[i].normal = outer_normal(dim, s.boundary[i], s.nodes, [xs[1] xs[2]])
        end
    else
        error("undef")
    end

    s.para["rho"] = rho
    s.para["mass"] = MathKits.get_volume(xs...) * rho
    s.para["I"] = I

    if dim == 1
        s.rotatable = false
    end
    return s
end

function outer_normal(dim::Int, c::Convex, nodes::Vector{Node}, xs::Array)
    normal = convex_normal(c, nodes, dim)
    nodesx = map(k->nodes[k].x0+nodes[k].d, c.link)
    xc = sum(nodesx)/length(nodesx)
    bias = 1.e-10
    xt = xc + normal*bias
    if pinpoly(xs, xt) == 1
        normal *= -1
    end
    return normal
end

function convex_normal(c, nodes, dim)
    if dim == 1
        normal = [1]
    elseif dim == 2
        normal = rotate_matrix(pi/2) * ((nodes[c.link[2]].x0+nodes[c.link[2]].d) - ((nodes[c.link[1]].x0+nodes[c.link[1]].d)))
    else
        error("undef")
    end
    return truncated_normalize(normal)
end

function truncated_normalize(v)
    v = normalize(v)
    a = map(x->x^2, v)
    p = sortperm(v)
    for i = 1:length(v)
        v[p[i]] = sqrt(sum(a) - sum(a[p[[(k==i ? false : true) for k=1:length(v)]]]))
    end
    return v
end

function set_fixed_u!(s::RigidStructure, u::Vector{Float64})
    s.u_is_fixed = true
    s.fixed_u = u
    for i = 1:s.nnp 
        s.nodes[i].u = copy(s.fixed_u)
        s.nodes[i].a = 0.
    end
    s.system.u = copy(s.fixed_u)
    s.system.a = zeros(Float64, s.dim)
end

function set_fixed_omega!(s::RigidStructure, omega::Float64)
    s.omega_is_fixed = true
    s.fixed_omega = omega
    s.system.omega = omega
    s.system.domega = 0.
end