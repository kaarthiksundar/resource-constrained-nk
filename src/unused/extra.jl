

function create_primal_model(data::Dict{String,Any}, model=Model())

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    ref = PMs.build_ref(data)[:nw][0]

    return model
end 

function create_dual_load_shedding_model(data::Dict{String,Any}, model=Model())

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    ref = PMs.build_ref(data)[:nw][0]

    # dual variables
    @variable(model, kcl[i in keys(ref[:bus])])
    @variable(model, pg_max[i in keys(ref[:bus])] >= 0)
    @variable(model, dc_lb[i in keys(ref[:branch])] >= 0)
    @variable(model, dc_ub[i in keys(ref[:branch])] >= 0)
    @variable(model, thermal_lb[i in keys(ref[:branch])] >= 0)
    @variable(model, thermal_ub[i in keys(ref[:branch])] >= 0)
    @variable(model, loadshed[i in keys(ref[:bus])] >= 0)
    @variable(model, ref_bus[i in keys(ref[:ref_buses])])

    # interdiction variables 
    @variable(model, x[i in keys(ref[:branch])], Bin)

    # auxiliary definitions for formulating the constraints
    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))

    # dual constraints 
    # (a) constraints corresponding to primal variable va
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        is_ref_bus = i in keys(ref[:ref_buses])
        if (is_ref_bus)
            @constraint(model, ref_bus[i] + sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dc_lb[l] - dc_ub[l]) for (l,f,t) in bus_arcs) == 0)
        else 
            @constraint(model, sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dc_lb[l] - dc_ub[l]) for (l,f,t) in bus_arcs) == 0)
        end
    end 

    # (b) constraints corresponding to primal variable pg
    for (i, gen) in ref[:gen]
        @constraint(model, - kcl[gen["gen_bus"]] - pg_max[i] == 0)
    end

    # (c) constraints corresponding to primal variable p
    for (l,i,j) in ref[:arcs_from]
        @constraint(model, kcl[i] - kcl[j] + thermal_lb[l] - thermal_ub[l] + dc_lb[l] - dc_ub[l] == 0)
    end

    # (d) constraints corresponding to primal variable ld
    for (i, bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        @constraint(model, - sum(load["pd"] for load in bus_loads) * kcl[i] - loadshed[i] <= sum(load["pd"] for load in bus_loads))
    end


    return model
end

function get_bounds(ref, model, variables)

    bounds = Dict{Symbol,Any}( [var => Dict{Any, Any}() for var in variables ])
    thermal_lb = getindex(model, :thermal_lb)
    thermal_ub = getindex(model, :thermal_ub)

    branches = keys(ref[:branch])

    for branch in branches
        @objective(model, Max, thermal_lb[branch])
        solve(model)
        bounds[:thermal_lb][branch] = getobjectivevalue(model)

        @objective(model, Max, thermal_ub[branch])
        solve(model)
        bounds[:thermal_ub][branch] = getobjectivevalue(model)
    end 

    return bounds
end 

function update_dual_bounds(p::Problem, c::Configuration)

    ref = get_ref(p)
    current_dual_bounds = get_dual_bounds(p)
    L = get_best_incumbent(p)
    U = get_upper_bound(p)

    m = Model()
    solver = GurobiSolver(LogToConsole=0)
    (use_cplex(c)) && (solver = CplexSolver())
    setsolver(m, solver)

    if get_problem_type(c) == :plain
        @variable(m, x[keys(ref[:branch])], Bin)
        @constraint(m, budget, sum(x) == get_k(c))
    elseif get_problem_type(c) == :topological 
        # connectivty enforced using single commodity flow
        @variable(m, x[i in keys(ref[:branch])], Bin)
        
        @variable(m, y[i in keys(ref[:bus])], Bin)
        @variable(m, f[(l,i,j) in ref[:arcs]] >= 0)
        @variable(m, dummy_x[i in keys(ref[:bus])], Bin)
        @variable(m, dummy_f[i in keys(ref[:bus])] >= 0)        
        
        @constraint(m, budget, sum(x) == get_k(c))
        @constraint(m, sum(dummy_x) == 1)
        
        # single commodity flow variable bounds
        @constraint(m, [(l,i,j) in ref[:arcs]], f[(l,i,j)] <= get_k(c) * x[l])
        @constraint(m, [i in keys(ref[:bus])], dummy_f[i] <= (get_k(c) + 1) * dummy_x[i])

        # connecting y and x variables
        @constraint(m, [i in keys(ref[:branch])], x[i] <= y[ref[:branch]["f_bus"]])
        @constraint(m, [i in keys(ref[:branch])], x[i] <= y[ref[:branch]["t_bus"]])

        # single commodity flow balance constraints
        @constraint(m, [k in keys(ref[:bus])], dummy_f[k] + sum(f[(l,j,i)] - f[(l,i,j)] for (l,i,j) in ref[:bus_arcs][k]) == y[k])
    else # planar
        spatial_map = get_spatial_map(p)
        @variable(m, x[keys(ref[branch])], Bin)
        @variable(m, y[i in keys(ref[:bus])], Bin)

        @constraint(m, budget, sum(x) <= get_k(c))
        @constraint(m, sum(y) == 1)
        @constraint(m, [i in keys(ref[:branch])], x[i] <= sum(y[j] for j in spatial_map[i]))
    end 

    # dual variables 
    @variable(m, kcl[i in keys(ref[:bus])])
    @variable(m, pg_max[i in keys(ref[:gen])] >= 0)
    @variable(m, dc_lb[i in keys(ref[:branch])] >= 0)
    @variable(m, dc_ub[i in keys(ref[:branch])] >= 0)
    @variable(m, thermal_lb[i in keys(ref[:branch])] >= 0)
    @variable(m, thermal_ub[i in keys(ref[:branch])] >= 0)
    @variable(m, loadshed[i in keys(ref[:load])] >= 0)
    @variable(m, ref_bus[i in keys(ref[:ref_buses])])

    # McCormick variables 
    @variable(m, thermal_lb_x[i in keys(ref[:branch])])
    @variable(m, thermal_ub_x[i in keys(ref[:branch])])

    # additional definitions for formulating the constraints
    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))
    ld_obj = Dict(i => 1.0 for i in keys(ref[:load]))
    for (i, load) in ref[:load]
        (load["pd"] < 0) && (ld_obj[i] = -1.0)
    end 

    # (a) constraints corresponding to primal variable va
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        is_ref_bus = i in keys(ref[:ref_buses])
        if (is_ref_bus)
            @constraint(m, ref_bus[i] + sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dc_lb[l] - dc_ub[l]) for (l,f,t) in bus_arcs) == 0)
        else 
            @constraint(m, sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dc_lb[l] - dc_ub[l]) for (l,f,t) in bus_arcs) == 0)
        end
    end 

    # (b) constraints corresponding to primal variable pg
    for (i, gen) in ref[:gen]
        @constraint(m, - kcl[gen["gen_bus"]] - pg_max[i] <= 0)
    end

    # (c) constraints corresponding to primal variable p
    for (l,i,j) in ref[:arcs_from]
        @constraint(m, kcl[i] - kcl[j] + thermal_lb[l] - thermal_ub[l] + dc_lb[l] - dc_ub[l] == 0)
    end

    # (d) constraints corresponding to primal variable ld
    for (i, load) in ref[:load]
        load_bus = load["load_bus"]
        @constraint(m, - load["pd"] * kcl[load_bus] - loadshed[i] <= load["pd"] * ld_obj[i])
    end 

    @expression(m, kcl_expr[i in keys(ref[:bus])], sum( -ref[:load][j]["pd"] for j in ref[:bus_loads][i] ) * kcl[i] )
    @expression(m, pgmax_expr, sum( -gen["pmax"] * pg_max[i] for (i,gen) in ref[:gen] ) )
    @expression(m, tmin_expr, sum( -branch["rate_a"] * thermal_lb[i] + branch["rate_a"] * thermal_lb_x[i] for (i,branch) in ref[:branch] ) )
    @expression(m, tmax_expr, sum( -branch["rate_a"] * thermal_ub[i] + branch["rate_a"] * thermal_ub_x[i] for (i,branch) in ref[:branch] ) )
    @expression(m, load_shed_expr, sum( -loadshed[i] for (i,load) in ref[:load] ) )


    @expression(m, obj_expr, sum(kcl_expr) + pgmax_expr + tmin_expr + tmax_expr + load_shed_expr)
    @constraint(m, L <= obj_expr)
    @constraint(m, obj_expr <= U)

    for (i,branch) in ref[:branch]
        add_reformulation(m, thermal_lb_x[i], x[i], thermal_lb[i], current_dual_bounds[i][:lb])
        add_reformulation(m, thermal_ub_x[i], x[i], thermal_ub[i], current_dual_bounds[i][:ub])
    end 
    
    for (i,branch) in ref[:branch]
        @objective(m, Max, thermal_lb[i])
        status = JuMP.solve(m)
        current_dual_bounds[i][:lb] = getobjectivevalue(m)
        xval = i -> getvalue(x[i])
        xvals = Dict( i => xval(i) for i in keys(p.ref[:branch]))
        selected_branches = collect( keys( filter( (i, val) -> val > 0.9, xvals ) ) )
        println(selected_branches)

        @objective(m, Max, thermal_ub[i])
        status = JuMP.solve(m)
        current_dual_bounds[i][:ub] = getobjectivevalue(m)
        xvals = Dict( i => xval(i) for i in keys(p.ref[:branch]))
        selected_branches = collect( keys( filter( (i, val) -> val > 0.9, xvals ) ) )
        println(selected_branches)
    end 
    println()
    for (i,branch) in ref[:branch]
        debug(logger, "$(current_dual_bounds[i])")
    end

    set_dual_bounds(p, current_dual_bounds)
    return 
end 


function add_reformulation(model, xy, x_bin, y_cont, y_ub) 
    @constraint(model, xy >= y_cont - (1-x_bin)*y_ub)
    @constraint(model, xy <= y_cont)
    @constraint(model, xy <= y_ub*x_bin)

    return
end


function get_dual_bounds(p::Problem, c::Configuration)

    bounds = Dict(
        "flow_bound" => Dict(), 
        "pf_bound" => Dict()
        )


    M = 1e5
    ref = get_ref(p)

    for (i, branch) in ref[:branch]
        bounds["flow_bound"][i] = Dict( :<= => NaN , :>= => NaN ) 
        bounds["pf_bound"][i] = Dict( :<= => NaN , :>= => NaN ) 
    end 

    m = Model()
    solver = GurobiSolver(LogToConsole=0)
    (use_cplex(c)) && (solver = CplexSolver())
    setsolver(m, solver)

    if get_problem_type(c) ==:plain
        @variable(m, x[keys(ref[:branch])], Bin)
        @constraint(m, budget, sum(x) == get_k(c))
    elseif get_problem_type(c) == :topological 
        # connectivty enforced using single commodity flow
        @variable(m, x[i in keys(ref[:branch])], Bin)
        
        @variable(m, y[i in keys(ref[:bus])], Bin)
        @variable(m, f[(l,i,j) in ref[:arcs]] >= 0)
        @variable(m, dummy_x[i in keys(ref[:bus])], Bin)
        @variable(m, dummy_f[i in keys(ref[:bus])] >= 0)        
        
        @constraint(m, budget, sum(x) == get_k(c))
        @constraint(m, sum(dummy_x) == 1)
        
        # single commodity flow variable bounds
        @constraint(m, [(l,i,j) in ref[:arcs]], f[(l,i,j)] <= get_k(c) * x[l])
        @constraint(m, [i in keys(ref[:bus])], dummy_f[i] <= (get_k(c) + 1) * dummy_x[i])

        # connecting y and x variables
        @constraint(m, [i in keys(ref[:branch])], x[i] <= y[ref[:branch]["f_bus"]])
        @constraint(m, [i in keys(ref[:branch])], x[i] <= y[ref[:branch]["t_bus"]])

        # single commodity flow balance constraints
        @constraint(m, [k in keys(ref[:bus])], dummy_f[k] + sum(f[(l,j,i)] - f[(l,i,j)] for (l,i,j) in ref[:bus_arcs][k]) == y[k])
    else # planar
        spatial_map = get_spatial_map(p)
        @variable(m, x[keys(ref[branch])], Bin)
        @variable(m, y[i in keys(ref[:bus])], Bin)

        @constraint(m, budget, sum(x) <= get_k(c))
        @constraint(m, sum(y) == 1)
        @constraint(m, [i in keys(ref[:branch])], x[i] <= sum(y[j] for j in spatial_map[i]))
    end 

    @variable(m, va[i in keys(ref[:bus])])
    
    @variable(m, 0 <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    
    @variable(m, p[(l, i, j) in ref[:arcs_from]])

    ld_max = Dict(i => 0.0 for i in keys(ref[:load]))
    ld_obj = Dict(i => 1.0 for i in keys(ref[:load]))
    for (i, load) in ref[:load]
        (load["pd"] > 0) && (ld_max[i] = 1.0)
        (load["pd"] < 0) && (ld_obj[i] = 0.0)
    end 
    @variable(m, 0 <= ld[i in keys(ref[:load])] <= ld_max[i])
    
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

    @objective(m, Max, sum(load["pd"] * ld_obj[i] * ld[i] for (i, load) in ref[:load]))

    time = @elapsed status = JuMP.solve(m, suppress_warnings=true)

    print_with_color(:red, "Starting dual value bounds MIP solve\n")
    println("Time taken     :  $(round(time; digits=2)) sec.")
    println("")   


    xval = i -> getvalue(x[i])
    xvals = Dict( i => xval(i) for i in keys(ref[:branch]))

    debug(logger, "interdiction plan:")
    for (i,value) in xvals 
        debug(logger, "x[$i] => $value")
    end
    
    print_with_color(:red, "Starting dual value bounds LP solve\n")
    
    bounds = get_dual_bounds(ref, bounds, xvals)


    print_with_color(:red, "Dual value bounds obtained\n") 

    return bounds
end 

function get_dual_bounds(ref::Dict, bounds::Dict, xvals::Dict)
    
    m = Model()
    M = 1e5
    solver = GurobiSolver(LogToConsole=0)
    setsolver(m, solver)

    constraint_refs = Dict( "flow_bound" => Dict() , "pf_bound" => Dict())
    for (i, branch) in ref[:branch]
        constraint_refs["flow_bound"][i] = Dict()
        constraint_refs["pf_bound"][i] = Dict()
    end 

    @variable(m, va[i in keys(ref[:bus])])
    
    @variable(m, 0 <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    
    @variable(m, p[(l, i, j) in ref[:arcs_from]])

    ld_max = Dict(i => 0.0 for i in keys(ref[:load]))
    ld_obj = Dict(i => 1.0 for i in keys(ref[:load]))
    for (i, load) in ref[:load]
        (load["pd"] > 0) && (ld_max[i] = 1.0)
        (load["pd"] < 0) && (ld_obj[i] = 0.0)
    end 
    @variable(m, 0 <= ld[i in keys(ref[:load])] <= ld_max[i])
    
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
        constraint_refs["pf_bound"][i][:>=] = @constraint(m, -M * xvals[i] <= p_fr + b*(va_fr - va_to))
        constraint_refs["pf_bound"][i][:<=] = @constraint(m, p_fr + b*(va_fr - va_to) <=  M * xvals[i])
    end

    for (l, i, j) in ref[:arcs_from]
        constraint_refs["flow_bound"][l][:>=] = @constraint(m, -ref[:branch][l]["rate_a"] * (1 - xvals[l]) <= p[(l, i, j)])
        constraint_refs["flow_bound"][l][:<=] = @constraint(m, p[(l, i, j)] <= ref[:branch][l]["rate_a"] * (1 - xvals[l]))
    end 

    @objective(m, Min, sum(load["pd"] * ld_obj[i] * ld[i] for (i, load) in ref[:load]))

    time = @elapsed status = JuMP.solve(m, suppress_warnings=true)

    for (i, branch) in ref[:branch]
        bounds["flow_bound"][i][:<=] = getdual(constraint_refs["flow_bound"][i][:<=])
        bounds["flow_bound"][i][:>=] = getdual(constraint_refs["flow_bound"][i][:>=])
        bounds["pf_bound"][i][:<=] = getdual(constraint_refs["pf_bound"][i][:<=])
        bounds["pf_bound"][i][:>=] = getdual(constraint_refs["pf_bound"][i][:>=])
    end

    return bounds
end 