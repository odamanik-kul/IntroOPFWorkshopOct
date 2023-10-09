function init_model!(m::Model, grid_data::Dict,time_steps::Int,load_time_series_file::String)

    ##### SETS
    m.ext[:sets] = Dict()
    ##### TIME STEPS
    m.ext[:sets][:time_steps] = collect(1:time_steps) # set of time steps

    # Buses
    M = m.ext[:sets][:bus] = [bus_id for (bus_id,bus) in grid_data["bus"]]
    MSL = m.ext[:sets][:bus_slack] = [bus_id for (bus_id,bus) in grid_data["bus"] if bus["bus_type"] == 3] # set of slack buses
    # Generators
    G = m.ext[:sets][:gen] = [gen_id for (gen_id,gen) in grid_data["gen"]]
    # Branch
    BR = m.ext[:sets][:branch] = [br_id for (br_id,br) in grid_data["branch"]] # branches
    # Topology
    arcs_fr = m.ext[:sets][:arcs_fr] = [(br_id,string(br["f_bus"]),string(br["t_bus"])) for (br_id,br) in grid_data["branch"]] # arc branch from to
    arcs_to = m.ext[:sets][:arcs_to] = [(br_id,string(br["t_bus"]),string(br["f_bus"])) for (br_id,br) in grid_data["branch"]] # arc branch to from
    arcs = m.ext[:sets][:arcs] = [arcs_fr;arcs_to]
    bus_arcs = Dict((b_id, Tuple{String,String,String}[]) for b_id in m.ext[:sets][:bus])
    for (l,i,j) in arcs
        push!(bus_arcs[i], (l,i,j))
    end
    m.ext[:sets][:bus_arcs] = bus_arcs
    
    # Create bus pairs
    m.ext[:sets][:bus_ij] = []
    m.ext[:sets][:bus_ji] = []
    contain_bus_ij = []
    for a in m.ext[:sets][:arcs_fr]
        if !((a[2],a[3]) in contain_bus_ij)
            push!(m.ext[:sets][:bus_ij], (a[2],a[3]))
            push!(m.ext[:sets][:bus_ji], (a[3],a[2]))     
            push!(contain_bus_ij, (a[2],a[3]))
        end
    end
    m.ext[:sets][:bus_ij_ji] = [m.ext[:sets][:bus_ij];m.ext[:sets][:bus_ji]]

    # Load
    D = m.ext[:sets][:load] = [load_id for (load_id,load) in grid_data["load"]] # set of loads

    ##### PARAMETERS
    m.ext[:parameters] = Dict()
    baseMVA = m.ext[:parameters][:baseMVA] = grid_data["baseMVA"]
   
    # Buses
    m.ext[:parameters][:bus] = Dict()
    m.ext[:parameters][:bus][:vm_min] = Dict(bus_id => bus["vmin"] for (bus_id,bus) in grid_data["bus"])
    m.ext[:parameters][:bus][:vm_max] = Dict(bus_id => bus["vmax"] for (bus_id,bus) in grid_data["bus"])
    m.ext[:parameters][:bus][:va_min] = Dict(bus_id => -pi for (bus_id,bus) in grid_data["bus"])
    m.ext[:parameters][:bus][:va_max] = Dict(bus_id => pi for (bus_id,bus) in grid_data["bus"])

    # Branches
    m.ext[:parameters][:branch] = Dict()
    m.ext[:parameters][:branch][:s_max_a] = Dict(br_id => br["rate_a"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:s_max_b] = Dict(br_id => br["rate_b"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:r] = Dict(br_id => br["br_r"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:g] = Dict(br_id => real(1/(br["br_r"] + br["br_x"]im)) for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:x] = Dict(br_id => br["br_x"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:b] = Dict(br_id => imag(1/(br["br_r"] + br["br_x"]im)) for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:g_sh_fr] = Dict(br_id => br["g_fr"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:g_sh_to] = Dict(br_id => br["g_to"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:b_sh_fr] = Dict(br_id => br["b_fr"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:b_sh_to] = Dict(br_id => br["b_to"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:ang_min] = Dict(br_id => br["angmin"] for (br_id,br) in grid_data["branch"])
    m.ext[:parameters][:branch][:ang_max] = Dict(br_id => br["angmax"] for (br_id,br) in grid_data["branch"])

    m.ext[:parameters][:bus][:ij_ang_max] = Dict()
    m.ext[:parameters][:bus][:ij_ang_min] = Dict()
    m.ext[:parameters][:bus][:ji_ang_max] = Dict()
    m.ext[:parameters][:bus][:ji_ang_min] = Dict()
    contain_bus_ij = []
    for a in m.ext[:sets][:arcs_fr]
        if !((a[2],a[3]) in contain_bus_ij)
            m.ext[:parameters][:bus][:ij_ang_max][(a[2],a[3])] = grid_data["branch"][a[1]]["angmax"]
            m.ext[:parameters][:bus][:ij_ang_min][(a[2],a[3])] = grid_data["branch"][a[1]]["angmin"]
            m.ext[:parameters][:bus][:ji_ang_max][(a[3],a[2])] = grid_data["branch"][a[1]]["angmax"]
            m.ext[:parameters][:bus][:ji_ang_min][(a[3],a[2])] = grid_data["branch"][a[1]]["angmin"]
            push!(contain_bus_ij, (a[2],a[3]))
        end
    end
    m.ext[:parameters][:bus][:ij_ji_ang_max] = merge!(m.ext[:parameters][:bus][:ij_ang_max],m.ext[:parameters][:bus][:ji_ang_max])
    m.ext[:parameters][:bus][:ij_ji_ang_min] = merge!(m.ext[:parameters][:bus][:ij_ang_min],m.ext[:parameters][:bus][:ji_ang_min])

    # Generators
    m.ext[:parameters][:gen] = Dict()
    m.ext[:parameters][:gen][:bus] =  Dict(gen_id => string(gen["gen_bus"]) for (gen_id,gen) in grid_data["gen"])
    m.ext[:parameters][:gen][:cost_linear] = Dict()
    for (gen_id,gen) in grid_data["gen"]
        if isempty(gen["cost"])
            m.ext[:parameters][:gen][:cost_linear][gen_id] = 0
        else
            m.ext[:parameters][:gen][:cost_linear][gen_id] = gen["cost"][gen["ncost"]-1]
        end
    end
    m.ext[:parameters][:gen][:p_max] = Dict(gen_id => gen["pmax"] for (gen_id,gen) in grid_data["gen"])
    m.ext[:parameters][:gen][:p_min] = Dict(gen_id => gen["pmin"] for (gen_id,gen) in grid_data["gen"])
    m.ext[:parameters][:gen][:q_max] = Dict(gen_id => gen["qmax"] for (gen_id,gen) in grid_data["gen"])
    m.ext[:parameters][:gen][:q_min] = Dict(gen_id => gen["qmin"] for (gen_id,gen) in grid_data["gen"])

    # Loads
    m.ext[:parameters][:load] = Dict()
    m.ext[:parameters][:load][:bus] = Dict(load_id => string(load["load_bus"]) for (load_id,load) in grid_data["load"])
    m.ext[:parameters][:load][:p_d] = Dict(load_id => load["pd"] for (load_id,load) in grid_data["load"])
    m.ext[:parameters][:load][:q_d] = Dict(load_id => load["qd"] for (load_id,load) in grid_data["load"])
    m.ext[:parameters][:load][:cost_curt] = Dict(load_id => 100000 for (load_id,load) in grid_data["load"])
    bus_load = m.ext[:parameters][:bus_load]= Dict(m_ => [] for m_ in m.ext[:sets][:bus])
    for (m_,load_) in bus_load
        for (load_id,load_bus) in m.ext[:parameters][:load][:bus]
            if load_bus==m_
                push!(load_,load_id)
            end
        end
        if isempty(load_)
            delete!(bus_load,m_)
        end
    end

    ##### TIME SERIES
    m.ext[:time_series] = Dict()
    # Loads
    load_pd_time_series_df = DataFrame(XLSX.readtable(load_time_series_file,"pd"))
    load_pd_time_series_dict_init = Dict(pairs(eachcol(load_pd_time_series_df)))
    load_qd_time_series_df = DataFrame(XLSX.readtable(load_time_series_file,"qd"))
    load_qd_time_series_dict_init = Dict(pairs(eachcol(load_qd_time_series_df)))
    m.ext[:time_series][:load] = Dict()
    m.ext[:time_series][:load][:p] = Dict(string(df_id) => df for (df_id,df) in load_pd_time_series_dict_init)
    m.ext[:time_series][:load][:q] = Dict(string(df_id) => df for (df_id,df) in load_qd_time_series_dict_init)

    return m
end