function fetch_data(s::RigidStructure, field)
    f = getfield(s.nodes[1], field)
    d = Array{eltype(f)}(undef, s.dim, s.nnp)
    for node in s.nodes 
        d[:,node.id] = getfield(node, field)
    end
    return d 
end

function copy_structure(s::RigidStructure)
    return RigidStructure(s.dim, s.nnp, s.ndof, copy_nodes(s.nodes), copy_system(s.system), copy_boundary(s.boundary), s.para, copy(s.ext_f), s.movable, s.rotatable, s.u_is_fixed, s.fixed_u, s.omega_is_fixed, s.fixed_omega)
end

function copy_nodes(nodes)
    nodes1 = Vector{Node}(undef, length(nodes))
    for i = 1:length(nodes)
        nodes1[i] = copy_node(nodes[i])
    end
    return nodes1
end

function copy_node(node)
    return Node(node.id, copy(node.x0), copy(node.d), copy(node.u), copy(node.a))
end

function copy_system(s::RigidSystem)
    return RigidSystem(copy(s.x0), copy(s.x), copy(s.d), copy(s.u), copy(s.a), s.theta, s.omega, s.domega)
end

function copy_boundary(b::Vector{Convex})
    b1 = Vector{Convex}(undef, length(b))
    for i = 1:length(b)
        b1[i] = copy_convex(b[i])
    end
    return b1
end

function copy_convex(c::Convex)
    return Convex(c.id, c.link, copy(c.normal))
end