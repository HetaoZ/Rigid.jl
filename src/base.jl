mutable struct Node
    id::Int
    x0::Vector{Float64}
    d::Vector{Float64}
    u::Vector{Float64}
    a::Vector{Float64}
end
Node(dim::Int) = Node(0, zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim))

mutable struct Convex
    id::Int # global ID, not the boundary's local ID
    link::Vector{Int}
    normal::Vector{Float64}
end
Convex() = Convex(0,Int[], Float64[])

mutable struct RigidSystem
    x0::Vector{Float64} # 中心点的初始位置。
    x::Vector{Float64} # 中心点，也是旋转中心，默认取结点平均值。
    d::Vector{Float64}
    u::Vector{Float64}
    a::Vector{Float64}
    theta::Float64 # 实时角度(默认初始角度为零)
    omega::Float64 # 角速度
    domega::Float64 # 角加速度    
end
RigidSystem(dim) = RigidSystem(zeros(Float64,dim),zeros(Float64,dim), zeros(Float64,dim), zeros(Float64,dim), zeros(Float64,dim), 0., 0., 0.)

mutable struct RigidStructure
    dim::Int
    nnp::Int
    ndof::Int
    nodes::Vector{Node}
    system::RigidSystem
    boundary::Vector{Convex}
    para::Dict
    ext_f::Vector{Float64}
    movable::Bool
    rotatable::Bool
    u_is_fixed::Bool
    fixed_u::Vector{Float64}
    omega_is_fixed::Bool
    fixed_omega::Float64
end
RigidStructure(dim, nnp) = RigidStructure(dim, nnp, dim*nnp, Vector{Node}(undef,nnp), RigidSystem(dim), Vector{Convex}(undef, nnp), Dict(), zeros(Float64,dim*nnp), true, true, false, zeros(Float64,dim), false, 0.)