
using JuMP, Gurobi, Plots

## PARAMETERS
m = 2;  # num. of linear vars
n = 1;  # num. of integer vars -- set to 1 for simplicity
num_constraints = 3; # num. of constraints
c = 6; # int obj
h = [7; 4]; # lin obj
A = [1 2 2]; # int constraints
G = [1 3; 0 2; 3 1]; # lin constraints
b = [12 9 10]; # RHS

## MODEL
model = Model(Gurobi.Optimizer)
@variable(model, x >= 0, Int)
@variable(model, y[1:m] >= 0)
@constraint(model, [k=1:num_constraints], A[k]*x + sum(G[k,j]*y[j] for j in 1:m) <= b[k])
@objective(model, Max, c*x + sum(h[j]*y[j] for j in 1:m))
optimize!(model);
model_rel = Model(Gurobi.Optimizer)    ## RELAXATION
@variable(model_rel, x_rel >= 0)
@variable(model_rel, y_rel[1:m] >= 0)
@constraint(model_rel, [k=1:num_constraints], A[k]*x_rel + sum(G[k,j]*y_rel[j] for j in 1:m) <= b[k])
@objective(model_rel, Max, c*x_rel + sum(h[j]*y_rel[j] for j in 1:m))
optimize!(model_rel);

## PRINT SOLUTIONS
print("opt_obj_MIP = "*string(objective_value(model))*"\n");
print("opt_obj_LPR = "*string(objective_value(model_rel))*"\n");
print("opt_sol_MIP = "*string([value.(x)], [value.(y)])*"\n");
print("opt_sol_LPR = "*string([value.(x_rel)], [value.(y_rel)])*"\n");
## ABSOLUTE GAP
print("abs_gap = "*string(objective_value(model_rel)-objective_value(model))*"\n");

eta = 1;
gap_1 = zeros(14);
for (ind_lambda, lambda) in enumerate(0:0.1:1.3)
    model_m = Model(Gurobi.Optimizer)
    @variable(model_m, x_m >= 0, Int)
    @variable(model_m, y_m[1:m] >= 0)
    @constraint(model_m, [k=1:num_constraints], A[k]*x_m + sum(G[k,j]*y_m[j] for j in 1:m) <= b[k] - eta*A[k] - lambda*G[k,1])
    @objective(model_m, Max, c*x_m + sum(h[j]*y_m[j] for j in 1:m))
    optimize!(model_m);
    model_l = Model(Gurobi.Optimizer)    ## RELAXATION
    @variable(model_l, x_l >= 0)
    @variable(model_l, y_l[1:m] >= 0)
    @constraint(model_l, [k=1:num_constraints], A[k]*x_l + sum(G[k,j]*y_l[j] for j in 1:m) <= b[k] - eta*A[k] - lambda*G[k,1])
    @objective(model_l, Max, c*x_l + sum(h[j]*y_l[j] for j in 1:m))
    optimize!(model_l);
    gap_1[ind_lambda] = objective_value(model_l) - objective_value(model_m);
end
gap_2 = zeros(14);
for (ind_lambda, lambda) in enumerate(0:0.1:1.3)
    model_m = Model(Gurobi.Optimizer)
    @variable(model_m, x_m >= 0, Int)
    @variable(model_m, y_m[1:m] >= 0)
    @constraint(model_m, [k=1:num_constraints], A[k]*x_m + sum(G[k,j]*y_m[j] for j in 1:m) <= b[k] - eta*A[k] - lambda*G[k,2])
    @objective(model_m, Max, c*x_m + sum(h[j]*y_m[j] for j in 1:m))
    optimize!(model_m);
    model_l = Model(Gurobi.Optimizer)    ## RELAXATION
    @variable(model_l, x_l >= 0)
    @variable(model_l, y_l[1:m] >= 0)
    @constraint(model_l, [k=1:num_constraints], A[k]*x_l + sum(G[k,j]*y_l[j] for j in 1:m) <= b[k] - eta*A[k] - lambda*G[k,2])
    @objective(model_l, Max, c*x_l + sum(h[j]*y_l[j] for j in 1:m))
    optimize!(model_l);
    gap_2[ind_lambda] = objective_value(model_l) - objective_value(model_m);
end

x_range = range(0, 1.3, length = 14);
plot(x_range, gap_1, xlabel="Value of λ", ylabel="Absolute Gap", label="Γ(b1)", 
    linecolor=:yellowgreen, markershape=:circle, markercolor=:yellowgreen, markerstrokecolor=:yellowgreen, 
        markerstrokewidth=:2, linewidth=2, markersize=4)
plot!(x_range, gap_2, label="Γ(b2)", linecolor=:purple, markershape=:xcross, markercolor=:purple, 
    markerstrokecolor=:purple, markerstrokewidth=:2, linewidth=2, markersize=2)
