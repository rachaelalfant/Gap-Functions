using JuMP
using Gurobi
    
let 
    # initialize 
    iter = 0;                                   # number of IP evaluations 
    b_MAX = 13;                                 # largest RHS
    RHS = collect(1:b_MAX);                     # vector of RHSs
    Gap = Array{Float64}(undef, length(RHS));   # gap between IP & LP for each RHS
    
    eta = 1;        # set eta
    c = [21; 11];   # set objective vector
    A = [7; 4];     # set constraint matrix
    n = length(A);  # number of columns of constraint matrix

    while true
        iter += 1; # update number of evaluations

        # columns for which periodicity applies changes at each iter
        index_list = Array{Int}(undef, 0); 
    
        # solve the IP with RHS b_MAX
        IP_model = Model(Gurobi.Optimizer) 
        @variable(IP_model, x[1:2]>=0, Int) #4
        @objective(IP_model, Max, sum(c[i]*x[i] for i in 1:n))
        @constraint(IP_model, sum(A[i]*x[i] for i in 1:n) <= b_MAX)
        optimize!(IP_model);
        zIP = objective_value(IP_model);
        optIP = value.(x);

        # solve the LP with RHS b_MAX
        LP_model = Model(Gurobi.Optimizer)
        @variable(LP_model, y[1:2]>=0) #4
        @objective(LP_model, Max, sum(c[i]*y[i] for i in 1:n))
        @constraint(LP_model, sum(A[i]*y[i] for i in 1:n) <= b_MAX)
        optimize!(LP_model);
        zLP = objective_value(LP_model);
        optLP = value.(y);

        Gap[b_MAX] = zLP - zIP; # gap at b_MAX

        # determine columns for which periodicity applies
        for j in 1:n
            if optIP[j] >= eta && optLP[j] >= eta
                push!(index_list, j);
            end
        end

        # RHSs for which periodicity applies
        if isempty([index_list]) == false
            alt_RHS = [b_MAX - A[j] for j in index_list];
            print("alt_RHS = "*string(alt_RHS)*"\n");
            alt_RHS = alt_RHS[alt_RHS .>0]; # filter out zero
            
            # update remaining RHSs    
            for num in alt_RHS
                Gap[num] = Gap[b_MAX];
                RHS = deleteat!(RHS,num); 
            end
        end
    
        pop!(RHS); # delete RHS we just solved for
    
        if isempty(RHS) == true # break out of loop when we have the gap at each RHS
            break
        end

        # update b_MAX
        b_MAX = last(RHS);
        
    end
    print("Gap = "*string(Gap)*"\n")            # output gap for each RHS
    print("Num. Evals. = "*string(iter)*"\n")   # output number of evaluations
end