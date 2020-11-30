module Rigid

# -------------------------------------------------
using LinearAlgebra, Statistics
using MathKits
const MK = MathKits
include("base.jl")
# -------------------------------------------------
const TEMPERAL_DIFF_SCHEME = "forward" # "central"

# -------------------------------------------------
mutable struct RigidMesh
    dim::Int
    n::Int
    neq::Int
    nodes::Vector{Node}
    IEN::Array{Int, 1}
    edge::Vector{IB}
    movable::Bool
    fixed_node_id::Int
end
RigidMesh(dim::Int, n::Int) = RigidMesh(dim, n, dim*n, Vector{Node}(undef, n), zeros(Int, n), IB[], true, 0)
RigidMesh(dim::Int) = RigidMesh(dim, 1)

mutable struct RigidStructure
    dim::Int
    mesh::RigidMesh
    ext::Externals
    rho::Float64 # 密度
    mass::Float64 # 整体质量
    I::Float64 # 转动惯量
    center0::Vector{Float64} # 中心点的初始位置。
    center::Vector{Float64} # 中心点，也是旋转中心，默认取结点平均值。
    d::Vector{Float64} # center位移
    u::Vector{Float64} # center速度
    a::Vector{Float64} # center加速度
    f::Vector{Float64} # center力
    theta0::Float64 # 初始角度
    theta::Float64 # 实时角度
    omega::Float64 # 角速度
    domega::Float64 # 角加速度
    rotatable::Bool
    movable::Bool
    u_is_fixed::Bool
    fixed_u::Vector{Float64}
end
RigidStructure(dim::Int, n::Int) = RigidStructure(dim, RigidMesh(dim, n), Externals(dim*n), 0., 0., 0., zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), 0., 0., 0., 0., true, true, false, zeros(Float64, dim))

"""
Generate a rigid structure.
"""
function RigidStructure(dim::Int, rho::Real, xs::Array{Float64}...; proto_shape::String = "rect3D")
    s = RigidStructure(dim, length(xs[1]))
    for i in eachindex(s.mesh.nodes)
        s.mesh.nodes[i] = Node(i, ones(Int,dim), [xs[k][i] for k in 1:dim], [xs[k][i] for k in 1:dim], zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), true)
    end
    s.mesh.IEN = [i for i in 1:length(xs[1])]
    s.rho = rho
    s.center0 = [mean(xs[k]) for k in 1:dim]
    s.center = copy(s.center0)
    if dim == 1
        s.mesh.edge = [IB([s.mesh.nodes[1]], [-1.0], false); IB([s.mesh.nodes[2]], [1.0], false)]
        s.mass = MK.get_volume(xs...) * rho
    elseif dim == 2
        s.mesh.edge = Array{IB}(undef, length(xs[1]))
        for i in eachindex(s.mesh.edge)
            s.mesh.edge[i] = IB(2)
            node1 = s.mesh.nodes[i]
            node2 = i==length(s.mesh.nodes) ? s.mesh.nodes[1] : s.mesh.nodes[i+1]
            s.mesh.edge[i].nodes = [node1, node2]
            s.mesh.edge[i].n = normalize(MK.rotate_matrix(pi/2) * (node2.x - node1.x))
            s.mesh.edge[i].ignored = false
        end
        @warn("You have to calculate the rotational inertia by yourself!")
        s.mass = MK.get_volume(xs...) * rho
    elseif dim == 3
        # 根据proto_shape决定处理方法
        error("undef dim")
    else
        @warn("You need to calculate the edge, mass, and rotational inertia by yourself!")
    end
    return s
end
# -------------------------------------------------
# 函数

function advance_solid!(s::RigidStructure, dt::Float64)
    #println("---- advance_solid ----")
    if s.movable
        apply_solver!(s, dt)
    end
end

function apply_solver!(s::RigidStructure, dt::Float64)
    solver_explicit!(s, dt)
end

function check_solid!(s::RigidStructure)
    check_edge!(s.mesh.edge)
    check_ext!(s.ext)
end

function Base.copy!(s::RigidStructure)
    s1 = RigidStructure(s.dim, s.mesh.n)
    s1.mesh = copy!(s.mesh)
    s1.mesh.edge = copy!(s.mesh.edge)
    s1.ext = copy!(s.ext)
    s1.rho = s.rho
    s1.mass = s.mass
    s1.I = s.I
    s1.center0 = copy(s.center0)
    s1.center = copy(s.center)
    s1.d = copy(s.d)
    s1.u = copy(s.u)
    s1.a = copy(s.a)
    s1.f = copy(s.f)
    s1.theta0 = s.theta0
    s1.theta = s.theta
    s1.omega = s.omega
    s1.domega = s.domega
    s1.rotatable = s.rotatable
    s1.movable = s.movable
    s1.u_is_fixed = s.u_is_fixed
    s1.fixed_u = copy(s.fixed_u)
    return s1
end

function Base.copy!(m::RigidMesh)
    m1 = RigidMesh(m.dim, m.n)
    m1.nodes= copy!(m.nodes)
    m1.IEN = copy(m.IEN)
    m1.edge = copy!(m.edge)
    m1.movable = m.movable
    m1.fixed_node_id = m.fixed_node_id
    return m1
end

function output!(frame::Int, time::Float64; filepath::String = "../outputdata/")
    if frame == 0
        open(filepath*"time.txt","w") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    else
        open(filepath*"time.txt","a") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    end
end

function output!(s::RigidStructure; varname::String = "mesh", frame::Int = 1, filepath::String = "outputdata/")
    if varname == "mesh"
        data = [node.x[k] for node in s.mesh.nodes, k in 1:s.dim]
        open(filepath*"solid_"*varname*"_"*string(frame+FRAME_BASE)*".txt","w") do file
            writedlm(file, data)
        end
    else
        error("undef varname")
    end
end

function outputfig!(frames::UnitRange; drawnode::Int = 0, drawcenter::Bool = true, filepath::String = "../../outputdata/")
    x = Float64[]
    for frame in frames
        data = readdlm(filepath*"solid_mesh_"*string(frame+FRAME_BASE)*".txt")
        if drawcenter
            push!(x, mean(data))
        else
            push!(x, data[drawnode])
        end
    end
    time = readdlm(filepath*"time.txt")[:, 2] 
    t = time[frames .+ 1]
    plot(t, x)
end

function output_node_history(k::Int; dir::Int = 1, filepath::String = "outputdata/", figpath::String = "outputfig/", ref::String = "/home/hetao/Projects/JuliaProj/CommonData/monasse2012_piston.txt")
    time = readdlm(filepath*"time.txt")[:,2]
    frame = readdlm(filepath*"time.txt")[:,1]
    data = Array{Float64}(undef, length(frame))
    i = 0
    for fr = frame
        i += 1
        dir_data = readdlm(filepath*"solid_mesh_"*string(Int(fr)+FRAME_BASE)*".txt")[:, dir]
        if k == 0
            data[i] = sum(dir_data) / length(dir_data)
        else
            data[i] = dir_data[k]
        end
    end
    figure(1, figsize=(8, 6))
    plot(time, data, color = "b")
    refdata = readdlm(ref)
    plot(refdata[:,1], refdata[:,2], color = "r")
    savefig(figpath*"node_history.png",dpi=100)
    close(1)
end

function set_fixed_speed!(s::RigidStructure, u::Vector{Float64})
    s.u_is_fixed = true
    s.fixed_u = u
    s.u = copy(s.fixed_u)
end

"""
刚体动力学。参见：https://zhuanlan.zhihu.com/p/45800345
"""
function solver_explicit!(s::RigidStructure, dt::Float64)
    s.f = [sum(s.ext.f[k:s.dim:end]) for k in 1:s.dim]
    for i in 1:s.mesh.n
        s.mesh.nodes[i].f = s.ext.f[s.dim*(i-1)+1:s.dim*i]
    end
    if s.rotatable && s.dim > 1
        if s.dim == 2
            M = 0.
            for node in s.mesh.nodes
                M += MK.crossproduct(node.x - s.center, node.f)
            end
            s.domega = M / s.I
            s.omega += s.domega * dt
            if TEMPERAL_DIFF_SCHEME == "central"
                s.theta += (s.omega - s.domega * dt + s.omega)/2 * dt
            elseif TEMPERAL_DIFF_SCHEME == "forward"
                s.theta += s.omega * dt
            else
                error("undef scheme")
            end
        else
            error("undef dim")
        end
    else
        s.domega = 0.
        s.omega = 0.
        s.theta = 0.
    end

    if s.u_is_fixed
        s.a = zeros(Float64, s.dim)
        s.u = s.fixed_u
    else
        # DEBUG_TEST
        # s.f += [50000.,]

        s.a = s.f ./ s.mass
        s.u += s.a * dt
    end
    
    # @warn "Restrictions of the lifted body"
    # if s.center[2] <= 0.05 
    #     s.u[2] = max(0., s.u[2])
    #     # s.a[2] = max(0., s.a[2])
    #     # s.f[2] = max(0., s.f[2])
    # end


    if TEMPERAL_DIFF_SCHEME == "central"
        s.d += (s.u - s.a*dt + s.u)/2 * dt
    elseif TEMPERAL_DIFF_SCHEME == "forward"
        s.d += s.u * dt
    else
        error("undef scheme")
    end
    
    s.center = s.center0 + s.d 
    
    if s.dim == 1
        for node in s.mesh.nodes
            node.a = s.a
            node.u = s.u
            node.d = s.d
            node.x = node.x0 + node.d
        end
    elseif s.dim == 2
        for node in s.mesh.nodes
            r = [MK.rotate_matrix(s.theta) * (node.x0 - s.center0); 0.] # 考虑旋转
            Omega = [0., 0., s.omega]
            dOmega = [0., 0., s.domega]
            node.u = s.u + MK.crossproduct(Omega, r)[1:2]
            node.a = s.a + MK.crossproduct(dOmega, r)[1:2] + MK.crossproduct(Omega, MK.crossproduct(Omega, r))[1:2]
            node.x = s.center + r[1:2]
            node.d = node.x - node.x0
        end
    else
        error("undef dim")
    end
    for b in s.mesh.edge
        for ii in eachindex(b.nodes)
            b.nodes[ii] = copy!(s.mesh.nodes[Tuple(b.nodes[ii].i)...])
        end
        update_normal!(b)
    end
    # println("-- advance_solid: 1 --")
    # println(s.mesh.edge[1].nodes[1].x)
end

function update_normal!(b::IB)
    if length(b.nodes) == 1
        # do nothing
    elseif length(b.nodes) == 2
        b.n = MK.get_normal(b.nodes[2].x - b.nodes[1].x)
    else
        error("undef dim")
    end
end

################################
end
