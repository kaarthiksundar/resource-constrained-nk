function populate_model(p::Problem, c::Configuration)
    
    # initializing solver
    solver = GurobiSolver(LogToConsole=0)
    (use_cplex(c)) && (solver = CplexSolver())

    setsolver(p.model, solver)

    if get_problem_type(c) == :plain
        @variable(p.model, x[keys(p.ref[:branch])], Bin)
        @constraint(p.model, budget, sum(x) == get_k(c))
        @variable(p.model, 0 <= eta <= 1e4)
        @objective(p.model, Max, eta)
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

        # eta and objective
        @variable(p.model, 0 <= eta <= 1e4)
        @objective(p.model, Max, eta)
    end 
    return
end 

function add_no_good_cut(p::Problem, c::Configuration)
    current_solution = get_current_solution(p)
    k = get_k(c)
    x = getindex(p.model, :x)
    branch_indexes = collect(keys(get_ref(p)[:branch]))
    n = length(collect(keys(get_ref(p)[:branch])))
    complement_current_solution = setdiff(branch_indexes, current_solution)
    @constraint(p.model, sum(x[i] for i in current_solution) + sum((1-x[i]) for i in complement_current_solution) <= n-1)
    @constraint(p.model, sum(x[i] for i in current_solution) <= k-1)
    
    return 
end 

function add_cutting_plane(p::Problem, sol::Dict{String,Any})
    current_solution = get_current_solution(p)
    x = getindex(p.model, :x)
    eta = getindex(p.model, :eta)
    load_shed = get_current_incumbent(p)
    va = sol["va"]
    flow = sol["branch_flow"]
    dual_bounds = get_dual_bounds(p)

    var = Any[] 
    coeff = Any[]
    constant = load_shed

    for (i, branch) in p.ref[:branch]
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        arc = (i, branch["f_bus"], branch["t_bus"])
        g, b = PMs.calc_branch_y(branch)
        if !(i in current_solution) 
            push!(var, x[i])
            flow_bound = (flow[arc] > 0) ? dual_bounds[i][:ub] : dual_bounds[i][:lb]
            flow_val = round(abs(flow[arc]) * flow_bound; digits=4)
            push!(coeff, flow_val)
        end 
    end 

    expr = AffExpr(var, coeff, constant)

    @constraint(p.model, eta <= expr)
    
    return 
end 