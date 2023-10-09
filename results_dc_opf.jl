# Get results
baseMVA = m_dc_opf.ext[:parameters][:baseMVA]
# Sets
set_G = m_dc_opf.ext[:sets][:gen]
set_M = m_dc_opf.ext[:sets][:bus]
set_arcs = m_dc_opf.ext[:sets][:arcs]
set_arcs_fr = m_dc_opf.ext[:sets][:arcs_fr]
set_arcs_to = m_dc_opf.ext[:sets][:arcs_to]
# Generators
gen_p = value.(m_dc_opf.ext[:variables][:gen_p])
# Buses
bus_ang = value.(m_dc_opf.ext[:variables][:bus_ang])
# AC branch power flows
p_b_ac = value.(m_dc_opf.ext[:variables][:p_b_ac])

# Plot generator active power
plot_markersize = 7
plot_fontfamily = "Computer Modern"
plot_titlefontsize = 20
plot_guidefontsize = 16
plot_tickfontsize = 12
plot_legendfontsize = 12
plot_size = (720,480)
plot_legend = false #:topright # :bottomright
plot_legend_column = 2
gen_p_all_mw = [gen_p[g,1] for g in set_G].*baseMVA
plot_gen_p = bar(gen_p_all_mw,
                    framestyle = :box,
                    legend= plot_legend,
                    # palette=cgrad(:default, length(gen_ids), categorical = true),
                    fontfamily=plot_fontfamily,
                    # background_color=:transparent,
                    foreground_color=:black,
                    titlefontsize = plot_titlefontsize,
                    guidefontsize = plot_guidefontsize,
                    tickfontsize = plot_tickfontsize,
                    legendfontsize = plot_legendfontsize,
                    size = plot_size,left_margin = (2,:mm),bottom_margin = (4,:mm))
xlabel!("Generator")
ylabel!("Active power [MW]")
xticks!([1:1:length(set_G);], set_G)
y_u_lim = maximum(gen_p_all_mw) + 50
ylims!(0,y_u_lim)
title!("Generator Dispatch")
Plots.svg(joinpath(@__DIR__,"results","dc_opf","plot_gen_p.svg"))

# Plot branch active power
plot_markersize = 7
plot_fontfamily = "Computer Modern"
plot_titlefontsize = 20
plot_guidefontsize = 16
plot_tickfontsize = 12
plot_legendfontsize = 12
plot_size = (720,480)
plot_legend = false #:topright # :bottomright
plot_legend_column = 2
p_b_ac_all_mw = [p_b_ac[(br,i,j),1] for (br,i,j) in set_arcs_fr].*baseMVA
plot_p_b_ac = bar(p_b_ac_all_mw,
                    framestyle = :box,
                    legend= plot_legend,
                    # palette=cgrad(:default, length(gen_ids), categorical = true),
                    fontfamily=plot_fontfamily,
                    # background_color=:transparent,
                    foreground_color=:black,
                    titlefontsize = plot_titlefontsize,
                    guidefontsize = plot_guidefontsize,
                    tickfontsize = plot_tickfontsize,
                    legendfontsize = plot_legendfontsize,
                    size = plot_size,left_margin = (2,:mm),bottom_margin = (4,:mm))
xlabel!("(branch,i,j)")
ylabel!("Active power [MW]")
xticks!([1:1:length(set_arcs_fr);], ["("*b*","*i*","*j*")" for (b,i,j) in set_arcs_fr])
title!("Branch flows")
Plots.svg(joinpath(@__DIR__,"results","dc_opf","plot_p_bij.svg"))






value.(m_dc_opf.ext[:variables][:gen_p])
value.(m_dc_opf.ext[:variables][:p_b_ac])
value.(m_dc_opf.ext[:variables][:load_p_curt])
value.(m_ac_opf.ext[:variables][:gen_p])
value.(m_ac_opf.ext[:variables][:gen_q])
load_p_curt = value.(m_ac_opf.ext[:variables][:load_p_curt])
load_q_curt = value.(m_ac_opf.ext[:variables][:load_q_curt])
load_p = m_ac_opf.ext[:time_series][:load][:p]
load_q = m_ac_opf.ext[:time_series][:load][:q]
load_p_rem = Dict(d => [load_p[d][t]-load_p_curt[d,t] for t in time_steps] for d in m_ac_opf.ext[:sets][:load])
load_q_rem = Dict(d => [load_q[d][t]-load_q_curt[d,t] for t in time_steps] for d in m_ac_opf.ext[:sets][:load])
load_pf = Dict(d => [m_ac_opf.ext[:time_series][:load][:p][d][t]./sqrt(m_ac_opf.ext[:time_series][:load][:p][d][1]^2+m_ac_opf.ext[:time_series][:load][:q][d][1]^2) for t in time_steps] for d in m_ac_opf.ext[:sets][:load])
load_pf_rem = Dict(d => [(load_p_rem[d][t])/sqrt((load_p_rem[d][t])^2+(load_q_rem[d][t])^2) for t in time_steps] for d in m_ac_opf.ext[:sets][:load])

value.(m_ac_opf.ext[:variables][:br_p])