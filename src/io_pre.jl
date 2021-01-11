"""
Generate a rigid structure in clockwise order.
"""
function create_rigid(dim::Int, rho::Real, I::Real, xs::Vector...)
    if dim == 1
        area = xs[1][2]-xs[1][1]
    elseif dim == 2
        area = MathKits.polygon_area(xs[1], xs[2])
    else
        error("undef dim")
    end
    if  area > 0.
        # do nothing
    elseif area < 0.
        xs = (reverse(xs[1]), reverse(xs[2]))
    else
        error("Area = 0")
    end

    s = RigidStructure(dim, length(xs[1]))

    for i in eachindex(s.nodes)
        s.nodes[i] = Node(i, [xs[k][i] for k in 1:dim], zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim))
    end
    
    s.system.x0 = [mean(xs[k]) for k in 1:dim]
    s.system.x = [mean(xs[k]) for k in 1:dim]
    s.system.theta = 0.
    s.system.omega = 0.
    s.system.domega = 0.

    create_boundary!(s)

    s.para["rho"] = rho
    s.para["mass"] = MathKits.get_volume(xs...) * rho
    s.para["I"] = I

    if dim == 1
        s.rotatable = false
    end
    return s
end

function create_boundary!(s)
    dim = s.dim
    if dim == 1
        s.boundary = [Convex(1,[1],[1.0]), Convex(2,[2],[1.0])]
        if s.nodes[1].x0[1] < s.nodes[2].x0[1]
            s.boundary[1].normal *= -1
        else
            s.boundary[2].normal *= -1
        end
    elseif dim == 2
        
        poly_coords = fetch_data(s, :x0)+fetch_data(s, :d)
        poly_xs = Matrix{Float64}(undef, s.nnp, 2)
        for i = 1:s.nnp
            poly_xs[i,:] = poly_coords[i]
        end

        for i in eachindex(s.boundary)
            s.boundary[i] = Convex(i, [i,i+1], [0.,0.])
        end
        s.boundary[end].link[2] = 1
        for i in eachindex(s.boundary)
            s.boundary[i].normal = outer_normal(dim, s.nodes[s.boundary[i].link], poly_xs)
        end
    else
        error("undef")
    end
end

function outer_normal(dim, nodes, xs)
    normal = - convex_normal(dim, nodes)

    nodesx = map(node->node.x0+node.d, nodes)
    xc = sum(nodesx)/length(nodesx)
    xt = xc + normal*1.e-10

    if dim == 1
        if MathKits.between(xt, nodes[1].x0+nodes[1].d, nodes[2].x0+nodes[2].d)
            normal *= -1
        end
    elseif dim == 2
        if pinpoly(xs, xt[1], xt[2]) == 1
            normal *= -1
        end

    else
        error("undef")
    end
    return normal
end

function convex_normal(dim, nodes)
    if dim == 1
        normal = [1]
    elseif dim == 2
        normal = rotate_matrix(pi/2) * ((nodes[2].x0+nodes[2].d) - ((nodes[1].x0+nodes[1].d)))
    else
        error("undef")
    end
    return normal/norm(normal)
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
        s.nodes[i].a = zeros(Float64, s.dim)
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