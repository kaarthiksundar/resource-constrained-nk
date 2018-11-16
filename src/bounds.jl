
function initialize_dual_bounds(p::Problem, c::Configuration, lines)

    M = 1
    (use_bt(c)) && (M = get_effective_load(p))
    ref = get_ref(p)
    bounds = Dict( [ i => Dict(:lb => M, :ub => M) for (i, branch) in ref[:branch] ] )

    for (l, branch) in ref[:branch]
        i = branch["f_bus"]
        j = branch["t_bus"]
        arc = (l, i, j)
        (arc in lines[:lb]) && (bounds[l][:lb] = 0.0)
        (arc in lines[:ub]) && (bounds[l][:ub] = 0.0)
    end 

    set_dual_bounds(p, bounds)
    return 
end 

