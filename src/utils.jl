function fetch_data(s::RigidStructure, field)
    f = getfield(s.nodes[1], field)
    d = Array{eltype(f)}(undef, s.dim, s.nnp)
    for node in s.nodes 
        d[:,node.id] = getfield(node, field)
    end
    return d 
end