
function solve_inner_problem(prob::Problem, c::Configuration, model=Model())

    # initializing solver
    solver = GurobiSolver(LogToConsole=0)
    (use_cplex(c)) && (solver = CplexSolver())

    setsolver(model, solver)

    # deactivate the branches in prob.data 
    for index in get_current_solution(prob)
        prob.data["branch"][string(index)]["br_status"] = 0
    end 

    # post the dc load shedding model - ignore the dc lines and theta diff bounds
    @assert !haskey(prob.data, "multinetwork")
    @assert !haskey(prob.data, "conductors")

    ref = PMs.build_ref(prob.data)[:nw][0]

    @variable(model, va[i in keys(ref[:bus])])
    
    @variable(model, 0 <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    
    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs_from]] <= ref[:branch][l]["rate_a"])

    ld_obj = Dict(i => 1.0 for i in keys(ref[:load]))
    for (i, load) in ref[:load]
        (load["pd"] < 0) && (ld_obj[i] = -1.0)
    end 
    @variable(model, 0 <= ld[i in keys(ref[:load])] <= 1.0)

    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))

    for (i, bus) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end 

    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]

        # Bus KCL
        @constraint(model, 
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
        @constraint(model, p_fr + b*(va_fr - va_to) == 0)
    end

    @objective(model, Min, sum(load["pd"] * ld_obj[i] * ld[i] for (i, load) in ref[:load]))

    # activate the branches in p.data 
    for index in get_current_solution(prob)
        prob.data["branch"][string(index)]["br_status"] = 1
    end 

    time = @elapsed status = JuMP.solve(model)

    for (l,i,j) in ref[:arcs_from]
        (abs(getdual(p[(l,i,j)])) > 1.0) && debug(logger, "($l,$i,$j) -> $(getvalue(p[(l,i,j)])), $(getdual(p[(l,i,j)]))")
    end 

    set_current_incumbent(prob, getobjectivevalue(model))

    if (get_current_incumbent(prob) - get_best_incumbent(prob) > 1e-4) 
        set_best_solution(prob, get_current_solution(prob))
        set_best_incumbent(prob, get_current_incumbent(prob))
        (get_problem_type(c) == :planar) && (set_center_bus_id(prob))
        for (i, bus) in ref[:bus]
            bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
            if (length(bus_loads) > 0)
                bus_load_shed = [getvalue(ld[load["index"]]) for load in bus_loads]
                bus_ld_obj = [ld_obj[load["index"]] for load in bus_loads]
                bus_pd = [load["pd"] for load in bus_loads]
                prob.bus_load_shed[i] = sum(bus_load_shed .* bus_ld_obj .* bus_pd)
            end
        end 
    end 

    update_opt_gap(prob)
    
    return time, status, model
end 

function get_inner_solution(prob::Problem, model)

    solution_dict = Dict{String,Any}()
    solution_dict["branch_flow"] = Dict{Any,Any}()
    solution_dict["va"] = Dict{Any,Any}()

    current_solution = get_current_solution(prob)
    p = getindex(model, :p)
    va = getindex(model, :va)

    ref = prob.ref
    
    for (l,i,j) in ref[:arcs_from]
        solution_dict["branch_flow"][(l,i,j)] = (l in current_solution) ? 0.0 : getvalue(p[(l,i,j)]) 
    end 

    for (i, bus) in ref[:bus]
        solution_dict["va"][i] = getvalue(va[i])
    end 

    return solution_dict
end 