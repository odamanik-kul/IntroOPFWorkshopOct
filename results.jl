

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