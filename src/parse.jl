using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--file", "-f"
        help = "pglib case file name"
        arg_type = String 
        default = "pglib_opf_case14_ieee.m"

        "--path", "-p"
        help = "data directory path"
        arg_type = String
        default = "../data/"

        "--type"
        help = "problem type [:plain, :topological, :planar]"
        arg_type = Symbol 
        default = :plain

        "--timeout", "-t"
        help = "time limit for the run in seconds"
        arg_type = Int 
        default = 86400

        "--gap", "-g"
        help = "optimality gap in %"
        arg_type = Float64
        default = 1e-2

        "-k" 
        help = "k value for interdiction"
        arg_type = Int 
        default = 2

        "--cplex"
        help = "flag for using cplex as the LP and MILP solver"
        action = :store_true

        "--parallel"
        help = "flag for performing parallel solves in dual variable bound computation"
        action = :store_true

        "--workers", "-w"
        help = "number of workers"
        arg_type = Int 
        default = 1

        "--log"
        help = "file logger flag"
        action = :store_true

        "--debug"
        help = "debug flag"
        action = :store_true

        "--heuristic"
        help = "run heuristic"
        action = :store_true 

        "--bt" 
        help = "bound tightening"
        action = :store_true
        
    end 

    return parse_args(s)
end