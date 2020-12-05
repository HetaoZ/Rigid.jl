function advance!(s::RigidStructure, dt::Float64, scheme::String)
    if s.movable
        if scheme == "explicit"
            explicit_solver!(s, dt)
        else
            error("undef")
        end
    end
end

"""
Ref：https://zhuanlan.zhihu.com/p/45800345
"""
function explicit_solver!(s::RigidStructure, dt::Float64)
    dim = s.dim
    f = [sum(s.ext_f[k:s.dim:end]) for k in 1:s.dim]
    
    
    if s.rotatable
        if s.omega_is_fixed
            s.system.omega = omega
            s.system.domega = 0.
        else
            if s.dim == 2
                M = 0.
                for node in s.nodes
                    M += cross_product(node.x0+node.d - s.system.x, s.ext_f[node.id*dim-dim+1:node.id*dim])[3]
                end
                s.system.domega = M / s.para["I"]
                s.system.omega += s.system.domega * dt
            else
                error("undef dim")
            end
        end
    else
        s.system.domega = 0.
        s.system.omega = 0.
    end

    if s.u_is_fixed
        s.system.u = copy(s.fixed_u)
        s.system.a = zeros(Float64, s.dim)
    else
        # DEBUG_TEST
        # s.f += [50000.,]

        s.system.a = f / s.para["mass"]
        s.system.u += s.system.a * dt
    end

    # @warn "Restrictions of the lifted body"
    # if s.center[2] <= 0.05 
    #     s.u[2] = max(0., s.u[2])
    #     # s.a[2] = max(0., s.a[2])
    #     # s.f[2] = max(0., s.f[2])
    # end

    if TEMPERAL_DIFF_SCHEME == "central"
        
    elseif TEMPERAL_DIFF_SCHEME == "forward"
        
    else
        error("undef scheme")
    end

    if TEMPERAL_DIFF_SCHEME == "central"
        s.system.theta += (s.system.omega - s.system.domega * dt + s.system.omega)/2 * dt
        s.system.d += (s.system.u - s.system.a*dt + s.system.u)/2 * dt
    elseif TEMPERAL_DIFF_SCHEME == "forward"
        s.system.theta += s.system.omega * dt
        s.system.d += s.system.u * dt
    else
        error("undef scheme")
    end
    
    s.system.x = s.system.x0 + s.system.d 
    
    if s.dim == 1
        for node in s.nodes
            node.a = s.system.a
            node.u = s.system.u
            node.d = s.system.d
        end
    elseif s.dim == 2
        for node in s.nodes
            r = [rotate_matrix(s.system.theta) * (node.x0 - s.system.x0); 0.] # 考虑旋转
            Omega = [0., 0., s.system.omega]
            dOmega = [0., 0., s.system.domega]
            node.u = s.system.u + cross_product(Omega, r)[1:2]
            node.a = s.system.a + cross_product(dOmega, r)[1:2] + cross_product(Omega, cross_product(Omega, r))[1:2]
            nodex = s.system.x + r[1:2]
            node.d = nodex - node.x0
        end
    else
        error("undef dim")
    end
    for i in 1:s.nnp
        s.boundary[i].normal = outer_normal(s.nodes[s.boundary[i].link], s.dim)
    end
end

function outer_normal(nodes, dim)
    normal = convex_normal(nodes, dim)
    nodesx = map(node->node.x0+node.d, nodes)
    xc = sum(nodesx)/length(nodesx)
    bias = 1.e-10
    xt = xc + normal*bias
    if pinpoly(xs, xt) == 1
        normal *= -1
    end
    return normal
end

function convex_normal(nodes, dim)
    if dim == 1
        normal = [1]
    elseif dim == 2
        normal = rotate_matrix(pi/2) * ((nodes[2].x0+nodes[2].d) - ((nodes[1].x0+nodes[1].d)))
    else
        error("undef")
    end
    return truncated_normalize(normal)
end