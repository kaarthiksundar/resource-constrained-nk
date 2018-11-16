
function run_heuristic(p::Problem, c::Configuration)

    ref = get_ref(p)
    load = ref[:load]
    gen = ref[:gen]
    bus_gens = ref[:bus_gens]
    bus_loads = ref[:bus_loads]
    gen_buses = [] 
    load_buses = [] 
    normal_buses = []
    bus_weights = Dict()

    for (i, bus) in ref[:bus]
        bus_load = (length(bus_loads[i]) > 0) ? sum([abs(load[j]["pd"]) for j in bus_loads[i]]) : 0.0
        bus_generation = (length(bus_gens[i]) > 0) ? sum([gen[j]["pmax"] for j in bus_gens[i]]) : 0.0
        effective_load = (bus_load > bus_generation) ? (bus_load - bus_generation) : 0.0
        if effective_load > 0.0
            push!(load_buses, i) 
            bus_weights[i] = effective_load
        elseif (effective_load == 0.0 && bus_generation > 0.0)
            push!(gen_buses, i)
            bus_weights[i] = 0.0
        else 
            push!(normal_buses, i)
            bus_weights[i] = 0.0
        end  
    end 

    m = Model() 
    solver = GurobiSolver(LogToConsole=0)
    solver = CplexSolver()
    (use_cplex(c)) && (solver = CplexSolver())

    setsolver(m, solver)

    if get_problem_type(c) ==:plain
        @variable(m, x[i in keys(ref[:branch])], Bin)
        @constraint(m, budget, sum(x) == get_k(c))
    elseif get_problem_type(c) == :topological 
        # connectivty enforced using single commodity flow
        @variable(m, x[i in keys(ref[:branch])], Bin)
        
        @variable(m, y1[i in keys(ref[:bus])], Bin)
        @variable(m, f[(l,i,j) in ref[:arcs]] >= 0)
        @variable(m, dummy_x[i in keys(ref[:bus])], Bin)
        @variable(m, dummy_f[i in keys(ref[:bus])] >= 0)        
        
        @constraint(m, budget, sum(x) == get_k(c))
        @constraint(m, sum(dummy_x) == 1)
        
        # single commodity flow variable bounds
        @constraint(m, [(l,i,j) in ref[:arcs]], f[(l,i,j)] <= get_k(c) * x[l])
        @constraint(m, [i in keys(ref[:bus])], dummy_f[i] <= (get_k(c) + 1) * dummy_x[i])

        # connecting y and x variables
        @constraint(m, [i in keys(ref[:branch])], x[i] <= y1[ref[:branch]["f_bus"]])
        @constraint(m, [i in keys(ref[:branch])], x[i] <= y1[ref[:branch]["t_bus"]])

        # single commodity flow balance constraints
        @constraint(m, [k in keys(ref[:bus])], dummy_f[k] + sum(f[(l,j,i)] - f[(l,i,j)] for (l,i,j) in ref[:bus_arcs][k]) == y1[k])

    else # planar
        @variable(m, x[i in keys(ref[:branch])], Bin)
        @constraint(m, budget, sum(x) == get_k(c))
    end 

    @variable(m, dummy_edges[i in gen_buses] == 0)
    @variable(m, y[i in keys(ref[:bus])], Bin)

    @constraint(m, gen_off[i in gen_buses], y[i] == 0)
    @objective(m, Max, sum(bus_weights[i]*y[i] for i in keys(ref[:bus])))


    bus_to_vertex_map = Dict() 
    vertex_to_bus_map = Dict()
    bus_to_vertex_map[0] = length(keys(ref[:bus])) + 1
    vertex_to_bus_map[length(keys(ref[:bus])) + 1] = 0
    count = 1
    for (i, bus) in ref[:bus]
        bus_to_vertex_map[i] = count 
        vertex_to_bus_map[count] = i 
        count += 1
    end 


    function lazycuts(cb)
        parallel_edges = Dict()
        x_vals = Dict( i => getvalue(x[i]) for i in keys(ref[:branch]) )
        y_vals = Dict( i => getvalue(y[i]) for i in load_buses )
        y_vals = filter!((i, val) -> val >= 1-1e-5, y_vals)
        g = SimpleWeightedGraph(length(keys(vertex_to_bus_map)))
        
        for (key, val) in x_vals
            f_bus = ref[:branch][key]["f_bus"]
            t_bus = ref[:branch][key]["t_bus"]
            i = bus_to_vertex_map[f_bus]
            j = bus_to_vertex_map[t_bus]
            (val == 0) && (val = 1e-5)
            add_edge!(g, i, j, val)
            if !haskey(parallel_edges, (i, j))
                parallel_edges[(i,j)] = Dict( :count => 1, :index => [key])
                parallel_edges[(j,i)] = Dict( :count => 1, :index => [key])
            else 
                parallel_edges[(i, j)][:count] += 1 
                push!(parallel_edges[(i,j)][:index], key)
                parallel_edges[(j,i)][:count] += 1
                push!(parallel_edges[(j,i)][:index], key)
            end 
        end 

        for gen in gen_buses
            f_bus = 0 
            t_bus = gen 
            i = bus_to_vertex_map[f_bus]
            j = bus_to_vertex_map[t_bus]
            add_edge!(g, i, j, 1e-5)
            if !haskey(parallel_edges, (i, j))
                parallel_edges[(i,j)] = Dict( :count => 1, :index => [0] )
                parallel_edges[(j,i)] = Dict( :count => 1, :index => [0] )
            end 
        end 

        r = dijkstra_shortest_paths(g, bus_to_vertex_map[0])
        for i in keys(y_vals)
            path = enumerate_paths(r, bus_to_vertex_map[i])
            @assert path[1] == length(keys(ref[:bus])) + 1 
            path_indexes = [] 
            for k in 2:length(path)-1 
                path_index = 0 
                index = parallel_edges[(path[k], path[k+1])]
                if index[:count] > 1 
                    path_index = (x_vals[index[:index][1]] > 1-1e-3) ? index[:index][2] : index[:index][1]
                else 
                    path_index = index[:index][1]
                end 
                push!(path_indexes, path_index)
            end

            if sum([x_vals[j] for j in path_indexes]) < y_vals[i]
                @lazyconstraint(cb, sum(x[j] for j in path_indexes) >= y[i])
            end 

        end 
    end 

    addlazycallback(m, lazycuts)
    status = JuMP.solve(m)

    x_vals = Dict(i => getvalue(x[i]) for i in keys(ref[:branch]))
    x_vals = filter((i, val) -> val > (1-1e-5), x_vals)
    y_vals = Dict( i => getvalue(y[i]) for i in keys(ref[:bus]) )
    y_vals = filter!((i, val) -> val >= 1-1e-5, y_vals)
    attack = Dict( :bus_ids => [], :branch_ids => [], :load_shed => 0.0 )
    for bus_id in collect(keys(y_vals))
        push!(attack[:bus_ids], bus_id)
    end
    for branch_id in collect(keys(x_vals))
        push!(attack[:branch_ids], branch_id)
    end
    attack[:load_shed] = getobjectivevalue(m)

    set_current_solution(p, attack[:branch_ids])
    return 
end 