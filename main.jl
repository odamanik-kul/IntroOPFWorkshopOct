using JuMP
using PowerModels
using Ipopt
using DataFrames
using XLSX
using CPLEX
const PM = PowerModels

# Setup file
case_file_name = "pglib_opf_case5_pjm.m"
case_file = joinpath(@__DIR__,"pg",case_file_name)
case_data = PM.parse_matpower(case_file)
load_time_series_file_name = "case5_pjm_load.xlsx"
load_time_series_file = joinpath(@__DIR__,"pg",load_time_series_file_name)

# General settings
time_steps = collect(1:1)

# Run DC power flow formulation
include("init_model.jl")
include("build_dc_opf.jl")
m_dc_opf = Model(Ipopt.Optimizer)
init_model!(m_dc_opf,case_data,1,load_time_series_file)
build_dc_opf!(m_dc_opf)
optimize!(m_dc_opf)

# Run AC power flow formulation using PowerModels
PM_dc_opf_results = PM.solve_dc_opf(case_file, Ipopt.Optimizer)

# Validation
println([objective_value(m_dc_opf),PM_dc_opf_results["objective"]])


# Run AC power flow formulation
include("init_model.jl")
include("build_ac_opf.jl")
m_ac_opf = Model(Ipopt.Optimizer)
init_model!(m_ac_opf,case_data,1,load_time_series_file)
build_ac_opf!(m_ac_opf)
optimize!(m_ac_opf)

# Run AC power flow formulation using PowerModels
PM_ac_opf_results = PM.solve_ac_opf(case_file, Ipopt.Optimizer)

# Validation
println([objective_value(m_ac_opf),PM_ac_opf_results["objective"]])

