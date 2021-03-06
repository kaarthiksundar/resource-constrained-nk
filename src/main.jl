using JuMP
using PowerModels
using Gurobi
using CPLEX
using Memento
using LightGraphs
using SimpleWeightedGraphs

if VERSION < v"0.7.0-"
    import Compat: round
end


PMs = PowerModels

Memento.setlevel!(getlogger("PowerModels"), "error")
logger = Memento.config!("info")

include("parse.jl")
include("utils.jl")
include("models.jl")
include("prob.jl")
include("solvers.jl")
include("bounds.jl")
include("heuristic.jl")
include("preprocess.jl")

# setting up the configuration for the run
configuration = Configuration(parse_commandline())
print_configuration(configuration)
(configuration.debug) && (Memento.setlevel!(logger, "debug"))

# populate problem data
problem = Problem(PMs.parse_file(get_case(configuration)))
if (get_problem_type(configuration) == :planar) 
    data_check(problem, configuration)
    populate_branch_lengths(problem)
    populate_spatial_map(problem, configuration)
end 

# preprocessing and heuristic runs
nontight_lines = Dict(:lb => [], :ub => [])
(use_bt(configuration)) && (nontight_lines = bt(problem, configuration))
(use_heuristic(configuration)) && (run_heuristic(problem, configuration))

set_effective_load(problem)
print_total_load(problem)
print_effective_load(problem)

# create and solve model
initialize_dual_bounds(problem, configuration, nontight_lines)
populate_model(problem, configuration)
table = Table() 
print_table_header(table) 
solve(problem, configuration, table)

write_result(problem, configuration)