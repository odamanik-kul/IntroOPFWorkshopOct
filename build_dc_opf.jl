function build_dc_opf!(m::Model)

    ##### TIME STEPS
    T = m.ext[:sets][:time_steps]
    
    ##### SETS
    # Buses
    M = m.ext[:sets][:bus] 
    MSL = m.ext[:sets][:bus_slack]
    # Generators
    G = m.ext[:sets][:gen]
    # Branch
    BR = m.ext[:sets][:branch]
    arcs_fr = m.ext[:sets][:arcs_fr]
    arcs_to = m.ext[:sets][:arcs_to]
    arcs = m.ext[:sets][:arcs]
    bus_arcs = m.ext[:sets][:bus_arcs]
    bus_ij_ji = m.ext[:sets][:bus_ij_ji]
    # Load
    D = m.ext[:sets][:load]

    ##### PARAMETERS
    baseMVA = m.ext[:parameters][:baseMVA]
    
    # Buses
    
    # Branches
    br_s_max_a = m.ext[:parameters][:branch][:s_max_a]
    br_s_max_b = m.ext[:parameters][:branch][:s_max_b]
    br_r = m.ext[:parameters][:branch][:r]
    br_g = m.ext[:parameters][:branch][:g]
    br_x = m.ext[:parameters][:branch][:x]
    br_b = m.ext[:parameters][:branch][:b]
    br_ang_min = m.ext[:parameters][:branch][:ang_min]
    br_ang_max = m.ext[:parameters][:branch][:ang_max]

    ij_ang_max = m.ext[:parameters][:bus][:ij_ang_max]
    ij_ang_max = m.ext[:parameters][:bus][:ij_ang_min]
    ij_ji_ang_max = m.ext[:parameters][:bus][:ij_ji_ang_max]
    ij_ji_ang_min = m.ext[:parameters][:bus][:ij_ji_ang_min]

    # Generators
    gen_bus = m.ext[:parameters][:gen][:bus]
    gen_ncost = m.ext[:parameters][:gen][:ncost]
    # Maximum ncost - To determine the cost function order
    max_gen_ncost = m.ext[:parameters][:gen][:max_ncost]
    gen_cost = m.ext[:parameters][:gen][:cost]
    gen_p_max = m.ext[:parameters][:gen][:p_max]
    gen_p_min = m.ext[:parameters][:gen][:p_min]
    gen_q_max = m.ext[:parameters][:gen][:q_max]
    gen_q_min = m.ext[:parameters][:gen][:q_min]

    # Loads
    load_bus = m.ext[:parameters][:load][:bus]
    load_p_d = m.ext[:parameters][:load][:p_d]
    load_q_d = m.ext[:parameters][:load][:q_d]
    bus_load = m.ext[:parameters][:bus_load]

    ##### TIME SERIES
    # Loads
    load_p = m.ext[:time_series][:load][:p]

    ##### VARIABLES
    m.ext[:variables] = Dict()
    # Generators
    gen_p = m.ext[:variables][:gen_p] = @variable(m, [G,T], base_name="gen_p")
    # Buses
    bus_ang = m.ext[:variables][:bus_ang] = @variable(m, [M,T], base_name="bus_ang")
    # AC branch power flows
    p_b_ac = m.ext[:variables][:p_b_ac] = @variable(m, [arcs,T], base_name="p_b_ac")
    
    ##### OBJECTIVE
    if max_gen_ncost == 1
        m.ext[:objective] = @objective(m, Min,
                sum(sum(gen_cost[g][1]
                        for g in G) for t in T)
        )
    elseif max_gen_ncost == 2
        m.ext[:objective] = @objective(m, Min,
                sum(sum(gen_cost[g][1]*gen_p[g,t] + gen_cost[g][2]
                        for g in G) for t in T)
        )
    elseif max_gen_ncost == 3
        m.ext[:objective] = @objective(m, Min,
                sum(sum(gen_cost[g][1]*gen_p[g,t]^2 + gen_cost[g][2]*gen_p[g,t] + gen_cost[g][3]
                        for g in G) for t in T)
        )
    elseif max_gen_ncost == 4
        m.ext[:objective] = @objective(m, Min,
                sum(sum(gen_cost[g][1]*gen_p[g,t]^3 + gen_cost[g][2]*gen_p[g,t]^2 + gen_cost[g][3]*gen_p[g,t] + + gen_cost[g][4]
                        for g in G) for t in T)
        )
    end

    ##### CONSTRAINTS
    m.ext[:constraints] = Dict()
    # Generator - Maximum and minimum active power limits
    m.ext[:constraints][:gen_p_max_lb] = @constraint(m, [g=G, t=T],
        gen_p_min[g] <= gen_p[g,t])
    m.ext[:constraints][:gen_p_max_ub] = @constraint(m, [g=G, t=T],
        gen_p[g,t] <= gen_p_max[g])

    # Active branch power flow limits
    m.ext[:constraints][:p_b_max_lb] = @constraint(m, [(br,i,j) = arcs, t=T],
        -br_s_max_a[br] <= p_b_ac[(br,i,j),t])
    m.ext[:constraints][:p_b_max_ub] = @constraint(m, [(br,i,j) = arcs, t=T],
        p_b_ac[(br,i,j),t] <= br_s_max_a[br])

    # Bus angle limits
    m.ext[:constraints][:ang_ij_lb] = @constraint(m, [(i,j) = bus_ij_ji, t=T],
        ij_ji_ang_min[(i,j)] <= bus_ang[i,t] - bus_ang[j,t])
    m.ext[:constraints][:ang_ij_ub] = @constraint(m, [(i,j) = bus_ij_ji, t=T],
        bus_ang[i,t] - bus_ang[j,t] <= ij_ji_ang_max[(i,j)])

    # Slack bus
    m.ext[:constraints][:ang_sl] = @constraint(m, [msl = MSL, t=T],
        bus_ang[msl,t] == 0)

    # Power flow constraints - AC power flow with DC approximation
    m.ext[:constraints][:p_b_ac] = @constraint(m, [(br,i,j) = arcs, t=T],
        p_b_ac[(br,i,j),t] == -br_b[br]*(bus_ang[i,t]-bus_ang[j,t]))

    # Nodal power balance - AC network
    m.ext[:constraints][:nodal_balance] = @constraint(m, [ma=M, t=T],
        sum(gen_p[gen_id,t] for (gen_id,gen) in gen_bus if gen == ma)
        - sum(load_p[load_id][t] for (load_id,load) in load_bus if load == ma)
        - sum(p_b_ac[(br,i,j),t] for (br,i,j) in bus_arcs[ma])
        == 0)

    return m
end