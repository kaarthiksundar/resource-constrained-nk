function bt(p::Problem, c::Configuration)

    ref = get_ref(p)

    m = Model()
    solver = GurobiSolver(LogToConsole=0)
    (use_cplex(c)) && (solver = CplexSolver())
    setsolver(m, solver)

    if get_problem_type(c) ==:plain
        @variable(m, x[i in keys(ref[:branch])], Bin)
        @constraint(m, budget, sum(x) == get_k(c))
    elseif get_problem_type(c) == :topological 
        @variable(p.model, x[keys(p.ref[:branch])], Bin)
        @variable(p.model, y[i in keys(p.ref[:bus])], Bin)
        @variable(p.model, f[(l,i,j) in p.ref[:arcs]] >= 0)
        @variable(p.model, dummy_x[i in keys(p.ref[:bus])], Bin)
        @variable(p.model, dummy_f[i in keys(p.ref[:bus])] >= 0)        
        
        @constraint(p.model, budget, sum(x) == get_k(c))
        @constraint(p.model, sum(dummy_x) == 1)
        
        # single commodity flow variable bounds
        @constraint(p.model, [(l,i,j) in p.ref[:arcs]], f[(l,i,j)] <= get_k(c) * x[l])
        @constraint(p.model, [i in keys(p.ref[:bus])], dummy_f[i] <= (get_k(c) + 1) * dummy_x[i])

        # connecting y and x variables
        @constraint(p.model, [i in keys(p.ref[:branch])], x[i] <= y[p.ref[:branch][i]["f_bus"]])
        @constraint(p.model, [i in keys(p.ref[:branch])], x[i] <= y[p.ref[:branch][i]["t_bus"]])

        # single commodity flow balance constraints
        @constraint(p.model, [k in keys(p.ref[:bus])], dummy_f[k] + sum(f[(l,j,i)] - f[(l,i,j)] for (l,i,j) in p.ref[:bus_arcs][k]) == y[k])
    else # planar
        @variable(m, x[i in keys(ref[:branch])], Bin)
        @constraint(m, budget, sum(x) == get_k(c))
    end 

    populate_bt_model(m, ref)
    pij = getindex(m, :p)
    bounds = Dict([(l, i, j) => Dict(:lb => NaN, :ub => NaN) for (l, i, j) in ref[:arcs_from]])
    non_tight_branches = Dict( :lb => [], :ub => [])

    for (l, i, j) in ref[:arcs_from]
        @objective(m, Min, pij[(l, i, j)]) 
        status = JuMP.solve(m, relaxation=true)
        bounds[(l, i, j)][:lb] = getobjectivevalue(m)
        (abs(bounds[(l, i, j)][:lb] + ref[:branch][l]["rate_a"]) > 1e-2) && (push!(non_tight_branches[:lb], (l, i, j)))
        print(*)

        @objective(m, Max, pij[(l,i,j)])
        status = JuMP.solve(m, relaxation=true)
        bounds[(l, i, j)][:ub] = getobjectivevalue(m)
        (abs(bounds[(l, i, j)][:ub] - ref[:branch][l]["rate_a"]) > 1e-2) && (push!(non_tight_branches[:ub], (l, i, j)))
        print(*)
    end 
    println()
    return non_tight_branches
end 

function populate_bt_model(m, ref)

    x = getindex(m, :x)
    M = 1e5

    @variable(m, va[i in keys(ref[:bus])])
    
    @variable(m, 0 <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    
    @variable(m, p[(l, i, j) in ref[:arcs_from]])

    ld_obj = Dict(i => 1.0 for i in keys(ref[:load]))
    for (i, load) in ref[:load]
        (load["pd"] < 0) && (ld_obj[i] = -1.0)
    end 
    @variable(m, 0 <= ld[i in keys(ref[:load])] <= 1.0)
    
    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))

    for (i, bus) in ref[:ref_buses]
        @constraint(m, va[i] == 0)
    end 

    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]

        # Bus KCL
        @constraint(m, 
            sum(p_expr[a] for a in ref[:bus_arcs][i]) == 
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] * (1-ld[load["index"]]) for load in bus_loads)
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PMs.calc_branch_y(branch)

        # Power flow constraints
        @constraint(m, -M * x[i] <= p_fr + b*(va_fr - va_to))
        @constraint(m, p_fr + b*(va_fr - va_to) <=  M * x[i])
    end

    for (l, i, j) in ref[:arcs_from]
        @constraint(m, -ref[:branch][l]["rate_a"] * (1 - x[l]) <= p[(l, i, j)])
        @constraint(m, p[(l, i, j)] <= ref[:branch][l]["rate_a"] * (1 - x[l]))
    end 

end 

