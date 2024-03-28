using JuMP
using Gurobi

# INITIALIZE PARAMETERS
iter = 0;                                   # number of IP evaluations 
bmax = 13;                                  # largest RHS
rhs = collect(1:bmax);                      # vector of RHSs
gap = Array{Float64}(undef, length(rhs));   # gap between IP & LP for each RHS

eta = 1;        # set eta
c = [21; 11];   # set objective vector
A = [7; 4];     # set constraint matrix
n = length(A);  # number of columns of constraint matrix

function periodic_gap(iter, bmax, rhs, gap, eta, c, A, n)
    while true
        iter += 1; # update number of evaluations

        # columns for which periodicity applies changes at each iter
        index_list = Array{Int}(undef, 0); 
    
        # solve the IP with RHS bmax
        IP_model = Model(Gurobi.Optimizer)
        set_silent(IP_model)
        @variable(IP_model, x[1:2]>=0, Int)
        @objective(IP_model, Max, sum(c[i]*x[i] for i in 1:n))
        @constraint(IP_model, sum(A[i]*x[i] for i in 1:n) <= bmax)
        optimize!(IP_model);
        zIP = objective_value(IP_model);
        optIP = value.(x);

        # solve the LP with RHS bmax
        LP_model = Model(Gurobi.Optimizer)
        set_silent(LP_model)
        @variable(LP_model, y[1:2]>=0)
        @objective(LP_model, Max, sum(c[i]*y[i] for i in 1:n))
        @constraint(LP_model, sum(A[i]*y[i] for i in 1:n) <= bmax)
        optimize!(LP_model);
        zLP = objective_value(LP_model);
        optLP = value.(y);

        gap[bmax] = zLP - zIP; # gap at bmax

        # determine columns for which periodicity applies
        for j in 1:n
            if optIP[j] >= eta && optLP[j] >= eta
                push!(index_list, j);
            end
        end

        # RHSs for which periodicity applies
        if isempty([index_list]) == false
            alt_rhs = [bmax - A[j] for j in index_list];
            print("alt_rhs = "*string(alt_rhs)*"\n");
            alt_rhs = alt_rhs[alt_rhs .>0]; # filter out zero
            
            # update remaining RHSs    
            for num in alt_rhs
                gap[num] = gap[bmax];
                rhs = deleteat!(rhs,num); 
            end
        end
    
        pop!(rhs); # delete RHS we just solved for

        if isempty(rhs) == true # break out of loop when we have the gap at each RHS
            break
        end

        # update bmax
        bmax = last(rhs);
        
    end
    return gap, iter;
end 

soln = periodic_gap(iter, bmax, rhs, gap, eta, c, A, n);
print("gap = "*string(soln[1])*"\n")            # output gap for each RHS
print("Num. Evals. = "*string(soln[2])*"\n")    # output number of evaluations
