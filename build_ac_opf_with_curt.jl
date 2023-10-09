function build_ac_opf!(m::Model)

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
    bus_vm_min = m.ext[:parameters][:bus][:vm_min]
    bus_vm_max = m.ext[:parameters][:bus][:vm_max]
    bus_va_min = m.ext[:parameters][:bus][:va_min]
    bus_va_max = m.ext[:parameters][:bus][:va_max]
   
    # Branches
    br_s_max_a = m.ext[:parameters][:branch][:s_max_a]
    br_s_max_b = m.ext[:parameters][:branch][:s_max_b]
    br_r = m.ext[:parameters][:branch][:r]
    br_g = m.ext[:parameters][:branch][:g]
    br_x = m.ext[:parameters][:branch][:x]
    br_b = m.ext[:parameters][:branch][:b]
    br_g_sh_fr = m.ext[:parameters][:branch][:g_sh_fr]
    br_g_sh_to = m.ext[:parameters][:branch][:g_sh_to]
    br_b_sh_fr = m.ext[:parameters][:branch][:b_sh_fr]
    br_b_sh_to = m.ext[:parameters][:branch][:b_sh_to]

    ij_ang_max = m.ext[:parameters][:bus][:ij_ang_max]
    ij_ang_max = m.ext[:parameters][:bus][:ij_ang_min]
    ij_ji_ang_max = m.ext[:parameters][:bus][:ij_ji_ang_max]
    ij_ji_ang_min = m.ext[:parameters][:bus][:ij_ji_ang_min]

    # Generators
    gen_bus = m.ext[:parameters][:gen][:bus]
    gen_cost_lin = m.ext[:parameters][:gen][:cost_linear]
    gen_p_max = m.ext[:parameters][:gen][:p_max]
    gen_p_min = m.ext[:parameters][:gen][:p_min]
    gen_q_max = m.ext[:parameters][:gen][:q_max]
    gen_q_min = m.ext[:parameters][:gen][:q_min]

    # Loads
    load_bus = m.ext[:parameters][:load][:bus]
    load_p_d = m.ext[:parameters][:load][:p_d]
    load_q_d = m.ext[:parameters][:load][:q_d]
    load_cost_curt = m.ext[:parameters][:load][:cost_curt]
    bus_load = m.ext[:parameters][:bus_load]

    ##### TIME SERIES
    # Loads
    load_p = m.ext[:time_series][:load][:p]
    load_q = m.ext[:time_series][:load][:q]

    ##### VARIABLES
    m.ext[:variables] = Dict()
    # Generators
    gen_p = m.ext[:variables][:gen_p] = @variable(m, [g=G,T], lower_bound=gen_p_min[g], upper_bound=gen_p_max[g], base_name="gen_p")
    gen_q = m.ext[:variables][:gen_q] = @variable(m, [g=G,T], lower_bound=gen_q_min[g], upper_bound=gen_q_max[g], base_name="gen_q")
    # Loads
    load_p_curt = m.ext[:variables][:load_p_curt] = @variable(m, [d=D,T], base_name="load_p_curt")
    load_q_curt = m.ext[:variables][:load_q_curt] = @variable(m, [d=D,T], base_name="load_q_curt")
    # Buses
    bus_v_a = m.ext[:variables][:bus_v_a] = @variable(m, [ma=M,T], lower_bound=bus_va_min[ma], upper_bound=bus_va_max[ma], base_name="bus_v_a")
    bus_v_m = m.ext[:variables][:bus_v_m] = @variable(m, [ma=M,T], lower_bound=bus_vm_min[ma], upper_bound=bus_vm_max[ma], base_name="bus_v_m")
    # AC branch power flows
    br_p = m.ext[:variables][:br_p] = @variable(m, [(br,i,j)=arcs,T], lower_bound=-br_s_max_a[br], upper_bound=br_s_max_a[br], base_name="br_p")
    br_q = m.ext[:variables][:br_q] = @variable(m, [(br,i,j)=arcs,T], lower_bound=-br_s_max_a[br], upper_bound=br_s_max_a[br], base_name="br_q")

    ##### OBJECTIVE FUNCTION
    m.ext[:objective] = @objective(m, Min,
            sum(sum(gen_cost_lin[g]*gen_p[g,t]
                    for g in G) for t in T)
            + sum(sum(load_cost_curt[d]*load_p_curt[d,t] for d in D) for t in T)
    )

    ##### CONSTRAINTS
    m.ext[:constraints] = Dict()

    # Branch power flow limits
    m.ext[:constraints][:br_s_max] = @NLconstraint(m, [(br,i,j)=arcs, t=T],
        br_p[(br,i,j),t]^2 + br_q[(br,i,j),t]^2 <= br_s_max_a[br]^2)

    # Bus angle difference limits
    m.ext[:constraints][:ang_ij_lb] = @constraint(m, [(i,j) = bus_ij_ji, t=T],
        ij_ji_ang_min[(i,j)] <= bus_v_a[i,t] - bus_v_a[j,t])
    m.ext[:constraints][:ang_ij_ub] = @constraint(m, [(i,j) = bus_ij_ji, t=T],
        bus_v_a[i,t] - bus_v_a[j,t] <= ij_ji_ang_max[(i,j)])

    # Slack bus
    m.ext[:constraints][:ang_sl] = @constraint(m, [msl = MSL, t=T],
        bus_v_a[msl,t] == 0)

    # Power flow constraints - AC power flow
    m.ext[:constraints][:p_bij] = @NLconstraint(m, [(br,i,j)=arcs_fr,t=T],
        br_p[(br,i,j),t] == 
        (br_g[br] + br_g_sh_fr[br])*bus_v_m[i,t]^2
        -(br_g[br] * bus_v_m[i,t] * bus_v_m[j,t] * cos(bus_v_a[i,t] - bus_v_a[j,t]))
        -(br_b[br] * bus_v_m[i,t] * bus_v_m[j,t] * sin(bus_v_a[i,t] - bus_v_a[j,t]))
    )
    m.ext[:constraints][:q_bij] = @NLconstraint(m, [(br,i,j)=arcs_fr,t=T],
        br_q[(br,i,j),t] ==
        -(br_b[br] + br_b_sh_fr[br])*bus_v_m[i,t]^2
        + (br_b[br] * bus_v_m[i,t] * bus_v_m[j,t] * cos(bus_v_a[i,t] - bus_v_a[j,t]))
        - (br_g[br] * bus_v_m[i,t] * bus_v_m[j,t] * sin(bus_v_a[i,t] - bus_v_a[j,t]))
    )
    m.ext[:constraints][:p_bji] = @NLconstraint(m, [(br,i,j)=arcs_to,t=T],
        br_p[(br,i,j),t] == 
        (br_g[br] + br_g_sh_to[br])*bus_v_m[i,t]^2
        -(br_g[br] * bus_v_m[i,t] * bus_v_m[j,t] * cos(bus_v_a[i,t] - bus_v_a[j,t]))
        -(br_b[br] * bus_v_m[i,t] * bus_v_m[j,t] * sin(bus_v_a[i,t] - bus_v_a[j,t]))
    )
    m.ext[:constraints][:q_bji] = @NLconstraint(m, [(br,i,j)=arcs_to,t=T],
        br_q[(br,i,j),t] ==
        -(br_b[br] + br_b_sh_to[br])*bus_v_m[i,t]^2
        + (br_b[br] * bus_v_m[i,t] * bus_v_m[j,t] * cos(bus_v_a[i,t] - bus_v_a[j,t]))
        - (br_g[br] * bus_v_m[i,t] * bus_v_m[j,t] * sin(bus_v_a[i,t] - bus_v_a[j,t]))
    )

    # Nodal power balance - Active power
    m.ext[:constraints][:nodal_p_balance] = @constraint(m, [ma=M, t=T],
        sum(gen_p[gen_id,t] for (gen_id,gen) in gen_bus if gen == ma)
        - sum(load_p[load_id][t] for (load_id,load) in load_bus if load == ma)
        + sum(load_p_curt[load_id,t] for (load_id,load) in load_bus if load == ma)
        - sum(br_p[(br,i,j),t] for (br,i,j) in bus_arcs[ma])
        == 0)
    # Nodal power balance - Reactive power
    m.ext[:constraints][:nodal_q_balance] = @constraint(m, [ma=M, t=T],
        sum(gen_q[gen_id,t] for (gen_id,gen) in gen_bus if gen == ma)
        - sum(load_q[load_id][t] for (load_id,load) in load_bus if load == ma)
        + sum(load_q_curt[load_id,t] for (load_id,load) in load_bus if load == ma)
        - sum(br_q[(br,i,j),t] for (br,i,j) in bus_arcs[ma])
        == 0)

    # Load curtailment limit
    m.ext[:constraints][:load_p_curt_ub] = @constraint(m, [d=D, t=T],
        load_p_curt[d,t] <= load_p[d][t])
    m.ext[:constraints][:load_p_curt_lb] = @constraint(m, [d=D, t=T],
        0 <= load_p_curt[d,t])
    m.ext[:constraints][:load_q_curt_ub] = @constraint(m, [d=D, t=T],
        load_q_curt[d,t] <= load_q[d][t])
    m.ext[:constraints][:load_q_curt_lb] = @constraint(m, [d=D, t=T],
        0 <= load_q_curt[d,t])
        
    m.ext[:constraints][:load_curt_pf] = @NLconstraint(m, [d=D, t=T],
        ((load_p[d][t]-load_p_curt[d,t])^2)*(load_p[d][t]^2 + load_q[d][t]^2) == (load_p[d][t]^2)*((load_p[d][t]-load_p_curt[d,t])^2+(load_q[d][t]-load_q_curt[d,t])^2)
    )

    return m
end