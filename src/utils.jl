
mutable struct Configuration
    config::Dict{Any,Any}
    case_filename::AbstractString
    case_path::AbstractString
    case_file::AbstractString
    problem_type::Symbol
    time_limit::Int
    opt_gap::Float64
    k::Int 
    d::Float64
    use_cplex::Bool 
    parallelize::Bool
    num_workers::Int 
    debug::Bool
    heuristic::Bool
    bt::Bool
    
    function Configuration(config)
        c = new()
        c.config = config 
        c.case_filename = config["file"]
        c.case_path = config["path"]
        c.case_file = c.case_path * c.case_filename
        c.problem_type = config["type"]
        c.time_limit = config["timeout"]
        c.opt_gap = config["gap"]
        c.k = config["k"]
        c.d = config["distance"]
        c.use_cplex = config["cplex"]
        c.parallelize = config["parallel"]
        c.num_workers = config["workers"]
        c.debug = config["debug"]
        c.heuristic = config["heuristic"]
        c.bt = config["bt"]
        return c
    end

end 

get_case(c::Configuration) = c.case_file 
get_problem_type(c::Configuration) = c.problem_type 
get_time_limit(c::Configuration) = c.time_limit
get_opt_gap(c::Configuration) = c.opt_gap
get_k(c::Configuration) = c.k 
get_d(c::Configuration) = c.d
use_cplex(c::Configuration) = c.use_cplex
use_heuristic(c::Configuration) = c.heuristic 
use_bt(c::Configuration) = c.bt

function print_configuration(c::Configuration)
    print_with_color(:red, "CLI arguments\n")
    println("Case           :  $(get_case(c))")
    println("Problem type   :  $(get_problem_type(c))")
    println("k value        :  $(get_k(c))")
    (get_problem_type(c) == :planar) && (println("d value        :  $(get_d(c)) km"))
    println("Time limit     :  $(get_time_limit(c))")
    println("Optimality gap :  $(get_opt_gap(c)*100)%")
    println("")
    return
end 

function get_distance(f_lat, f_lon, t_lat, t_lon)
    R = 6371e3 # radius of earth in meters
    phi_1 = deg2rad(f_lat)
    phi_2 = deg2rad(t_lat)
    del_phi = deg2rad(t_lat - f_lat)
    del_lambda = deg2rad(t_lon - f_lon)

    a = sin(del_phi/2) * sin(del_phi/2) + cos(phi_1) * cos(phi_2) * sin(del_lambda/2) * sin(del_lambda/2)
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    d = R * c 

    return d/1000
end

function sample_branch(f_lat, f_lon, t_lat, t_lon, d)
    fractions = collect(linspace(0.0, 1.0, 100))
    samples = []
    for f in fractions 
        (f == 0.0) && (push!(samples, (f_lat, f_lon)); continue)
        (f == 1.0) && (push!(samples, (t_lat, t_lon)); continue)
        A = sin((1-f)*d)/sin(d)
        B = sin(f*d)/sin(d)
        x = A * cos(deg2rad(f_lat)) * cos(deg2rad(f_lon)) + B * cos(deg2rad(t_lat)) * cos(deg2rad(t_lon))
        y = A * cos(deg2rad(f_lat)) * sin(deg2rad(f_lon)) + B * cos(deg2rad(t_lat)) * sin(deg2rad(t_lon))
        z = A * sin(deg2rad(f_lat)) + B * sin(deg2rad(t_lat))
        lat = atan2(z, sqrt(x*x+y*y))
        lon = atan2(y, x)
        push!(samples, (rad2deg(lat), rad2deg(lon)))
    end     
    return samples
end 

mutable struct Problem 
    data::Dict{String,Any}
    ref::Dict{Any,Any}
    spatial_map::Dict{Any,Any}
    total_load::Float64
    effective_load::Float64
    model::JuMP.Model 
    best_solution::Vector{Int}
    current_solution::Vector{Int}
    upper_bound::Float64
    best_incumbent::Float64
    current_incumbent::Float64
    opt_gap::Float64
    iterations::Int
    dual_bounds::Dict
    bus_ids::Vector{Int}
    
    function Problem(data)
        @assert !haskey(data, "multinetwork")
        @assert !haskey(data, "conductors")
        p = new()
        p.data = data
        p.ref = PMs.build_ref(data)[:nw][0]
        p.spatial_map = Dict()
        p.total_load = sum([abs(load["pd"]) for (i, load) in p.ref[:load]]) 
        p.effective_load = 0.0
        p.model = Model()
        p.best_solution = Vector{Int}()
        p.current_solution = Vector{Int}()
        p.upper_bound = Inf
        p.best_incumbent = -Inf
        p.current_incumbent = -Inf
        p.opt_gap = Inf
        p.iterations = 0
        p.dual_bounds = Dict()
        p.bus_ids = Vector{Int}()
        return p
    end 
end

function populate_branch_lengths(p::Problem)
    for (i, branch) in p.ref[:branch] 
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_lat = p.ref[:bus][f_bus]["latitude"]
        f_lon = p.ref[:bus][f_bus]["longitude"]
        t_lat = p.ref[:bus][t_bus]["latitude"]
        t_lon = p.ref[:bus][t_bus]["longitude"]
        distance = get_distance(f_lat, f_lon, t_lat, t_lon)
        p.ref[:branch][i]["length"] = distance
    end 
    return
end 

function populate_spatial_map(p::Problem, c::Configuration)
    ref = get_ref(p)
    spatial_map = Dict([i => [] for i in keys(ref[:branch])])

    for (i, branch) in ref[:branch]
        D = get_d(c)
        d = branch["length"]
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        f_lat = p.ref[:bus][f_bus]["latitude"]
        f_lon = p.ref[:bus][f_bus]["longitude"]
        t_lat = p.ref[:bus][t_bus]["latitude"]
        t_lon = p.ref[:bus][t_bus]["longitude"]
        sampled_branch = sample_branch(f_lat, f_lon, t_lat, t_lon, d)
        for (j, bus) in ref[:bus]
            bus_lat = bus["latitude"]
            bus_lon = bus["longitude"]
            for (lat, lon) in sampled_branch
                if (get_distance(bus_lat, bus_lon, lat, lon) <= get_d(c)/2)
                    push!(spatial_map[i], j)
                    break
                end 
            end 
        end 
    end 
    set_spatial_map(p, spatial_map)
    return 
end

function print_total_load(p::Problem)
    println("Total load     :  $(round(get_total_load(p); digits=2))")
    println("")
    return
end 

function print_effective_load(p::Problem)
    println("Effective load :  $(round(get_effective_load(p); digits=2))")
    println("")
    return
end 

get_data(p::Problem) = p.data
get_ref(p::Problem) = p.ref 
get_spatial_map(p::Problem) = p.spatial_map
get_total_load(p::Problem) = p.total_load
get_effective_load(p::Problem) = p.effective_load
get_model(p::Problem) = p.model 
get_best_solution(p::Problem) = p.best_solution 
get_current_solution(p::Problem) = p.current_solution
get_upper_bound(p::Problem) = p.upper_bound 
get_best_incumbent(p::Problem) = p.best_incumbent
get_current_incumbent(p::Problem) = p.current_incumbent
get_opt_gap(p::Problem) = p.opt_gap 
get_iteration_count(p::Problem) = p.iterations
get_dual_bounds(p::Problem) = p.dual_bounds
get_bus_ids(p::Problem) = p.bus_ids

function set_spatial_map(p::Problem, spatial_map)
    p.spatial_map = spatial_map
    return 
end 

function set_upper_bound(p::Problem, ub) 
    p.upper_bound = ub
    return 
end

function set_current_incumbent(p::Problem, lb) 
    p.current_incumbent = lb
    return 
end

function set_best_incumbent(p::Problem, lb) 
    p.best_incumbent = lb
    return 
end

function update_opt_gap(p::Problem)
    p.opt_gap = (p.upper_bound-p.best_incumbent)/(p.upper_bound+1e-5)
    return
end 

function increment_iteration_count(p::Problem)
    p.iterations += 1
    return
end

function set_current_solution(p::Problem)
    x = getindex(p.model, :x)
    xval = i -> getvalue(x[i])
    xvals = Dict( i => xval(i) for i in keys(p.ref[:branch]))
    selected_branches = collect( keys( filter( (i, val) -> val > 0.9, xvals ) ) )
    p.current_solution = selected_branches
    return 
end 

function set_current_solution(p::Problem, solution)
    p.current_solution = solution
    return 
end 

function set_best_solution(p::Problem, solution::Vector{Int})
    p.best_solution = solution
    return 
end 

function set_dual_bounds(p::Problem, dual_bounds::Dict)
    p.dual_bounds = dual_bounds
    return 
end

function set_effective_load(p::Problem)
    ref = get_ref(p)
    bus_gens = ref[:bus_gens]
    bus_loads = ref[:bus_loads]
    gens = ref[:gen]
    loads = ref[:load]
    total_load = get_total_load(p)
    effective_load = 0.0
    for (i, bus) in ref[:bus]
        bus_generation = (length(bus_gens[i]) > 0) ? sum([gens[j]["pmax"] for j in bus_gens[i]]) : 0.0
        bus_load = (length(bus_loads[i]) > 0) ? sum([abs(loads[j]["pd"]) for j in bus_loads[i]]) : 0.0
        (bus_load > bus_generation) && (effective_load += (bus_load - bus_generation))
    end 
    p.effective_load = effective_load
    return 
end 

function data_check(p::Problem, c::Configuration) 
    if get_problem_type(c) == :planar
        for (i, bus) in get_ref(p)[:bus]
            !(haskey(bus, "latitude")) && (error("lat-lon information missing and problem type specified as :planar"))
            !(haskey(bus, "longitude")) && (error("lat-lon information missing and problem type specified as :planar"))
        end
    end 
    return 
end 

mutable struct Result 
    status::Symbol 
    objective::Float64 
    solution::Vector{Int}
    bound::Float64
    time::Float64 
    opt_gap::Float64
    bus_ids::Vector{Int}
    
    function Result() 
        r = new()
        r.status = :Unknown
        r.objective = NaN 
        r.solution = Vector{Int}()
        r.bound = NaN 
        r.time = NaN 
        r.opt_gap = NaN
        r.bus_ids = Vector{Int}()
        return r
    end 
end 

get_status(r::Result) = r.status 
get_objective(r::Result) = r.objective
get_solution(r::Result) = r.solution 
get_bound(r::Result) = r.bound 
get_time(r::Result) = r.time 
get_opt_gap(r::Result) = r.opt_gap 
get_bus_ids(r::Result) = r.bus_ids

function set_status(r::Result, status::Symbol)
    r.status = status 
    return 
end 

function set_objective(r::Result, obj::Float64)
    r.objective = obj 
    return 
end 

function set_bound(r::Result, bound::Float64)
    r.bound = bound 
    return 
end 

function set_solution(r::Result, sol::Vector{Int})
    r.solution = sol 
    return 
end 

function set_time(r::Result, time::Float64) 
    r.time = time 
    return 
end 

function set_opt_gap(r::Result, gap::Float64)
    r.opt_gap = gap
    return 
end 

function set_bus_ids(r::Result, bus_ids::Vector{Int})
    r.bus_ids = bus_ids
    return 
end 

mutable struct Table
    fields::Dict{Symbol,Any}
    field_chars::Dict{Symbol,Any}
    all_fields::Vector{Symbol}
    total_field_chars::Int

    function Table()
        t = new()
        t.fields = Dict{Symbol,Any}()
        t.field_chars = Dict{Symbol,Any}()
        t.fields[:Incumbent] = NaN 
        t.fields[:BestBound] = NaN
        t.fields[:Gap] = NaN 
        t.fields[:Time] = 0.0
        t.fields[:Iter] = 0
        t.field_chars[:Incumbent] = 28
        t.field_chars[:BestBound] = 28
        t.field_chars[:Gap] = 12
        t.field_chars[:Time] = 12
        t.field_chars[:Iter] = 12
        t.all_fields = [:Iter, :Incumbent, :BestBound, :Gap, :Time]
        t.total_field_chars = 28*2 + 12*3
        return t
    end 
end 

function update_table(p::Problem, t::Table, time::Float64)
    t.fields[:Incumbent] = get_best_incumbent(p)
    t.fields[:BestBound] = get_upper_bound(p)
    t.fields[:Gap] = get_opt_gap(p)
    t.fields[:Time] += time 
    t.fields[:Iter] = get_iteration_count(p)
    return 
end

function print_table_header(t::Table)
    println("")
    println(repeat("-", t.total_field_chars))
    println(get_table_header(t))
    println(repeat("-", t.total_field_chars))
    return 
end

function print_table_line(t::Table, c::Configuration) 
    println(get_table_line(t, c))
    return 
end 

function print_table_footer(t::Table)
    println(repeat("-", t.total_field_chars))
    return 
end

function get_table_header(t::Table)
    line = ""
    for f in t.all_fields
        name = (f == :BestBound) ? "Best Bound" : string(f)
        (f == :Gap) && (name = "Gap (%)")
        padding = t.field_chars[f] - length(name)
        line *= repeat(" ", trunc(Int, floor(padding/2)))
        line *= name 
        line *= repeat(" ", trunc(Int, ceil(padding/2)))
    end 
    return line
end 

function get_table_line(t::Table, c::Configuration)
    line = ""
    for f in t.all_fields
        if f == :Iter 
            value = string(t.fields[f]) 
            padding = t.field_chars[f] - length(value)
            line *= repeat(" ", trunc(Int, floor(padding/2)))
            line *= value
            line *= repeat(" ", trunc(Int, ceil(padding/2)))
        end 

        if f == :Incumbent || f == :BestBound
            value = isinf(t.fields[f]) ? "-" : string(round(t.fields[f]; digits=2))
            padding = t.field_chars[f] - length(value)
            line *= repeat(" ", trunc(Int, floor(padding/2)))
            line *= value
            line *= repeat(" ", trunc(Int, ceil(padding/2)))
        end 

        if f == :Gap
            value = ""
            if isnan(t.fields[f])
                value = "-"
            elseif  isinf(t.fields[f])
                value = "∞"
            else
                value = round(t.fields[f]*100; digits=1)
                if length(string(value)) < t.field_chars[f]
                    if value < get_opt_gap(c)*100
                        value = "opt"
                    elseif value > 1000
                        value = "∞"
                    else
                        value = string(value)
                    end
                else
                    value = string(value)
                end
            end
            padding = t.field_chars[f] - length(value)
            line *= repeat(" ", trunc(Int, floor(padding/2)))
            line *= value
            line *= repeat(" ", trunc(Int, ceil(padding/2)))
        end 

        if f == :Time 
            value = string(round(t.fields[f]; digits=1))
            padding = t.field_chars[f] - length(value)
            line *= repeat(" ", trunc(Int, floor(padding/2)))
            line *= value
            line *= repeat(" ", trunc(Int, ceil(padding/2)))
        end 
    end 

    return line 
end 