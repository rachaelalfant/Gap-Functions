using JuMP, Gurobi

# INITIALIZE PARAMETERS
m = 11;
n = 4;
param_UB = 50;
param_LB = 30;
f = rand(param_LB:param_UB,n);
c = rand(1:param_UB,n,m);
d = rand(1:param_UB,m);

function facility_loc(m, n, param_UB, param_LB, f, c, d)
    ## WEAK MODEL
    model_w = Model(Gurobi.Optimizer)
    @variable(model_w, x_w[1:n], Bin)
    @variable(model_w, z_w[1:n,1:m], Bin)
    @constraint(model_w, [j=1:m], sum(z_w[i,j] for i in 1:n) == 1)
    @constraint(model_w, [i=1:n], sum(z_w[i,j] for j in 1:m) <= m*x_w[i])
    @objective(model_w, Min, sum(f[i]*x_w[i] for i in 1:n) + sum((c[i,j]*d[j]*z_w[i,j] for j in 1:m) for i in 1:n))
    optimize!(model_w);
    model_w_LP = Model(Gurobi.Optimizer)    ## RELAXATION
    @variable(model_w_LP, 0 <= x_w_LP[1:n] <= 1)
    @variable(model_w_LP, 0 <= z_w_LP[1:n,1:m] <= 1)
    @constraint(model_w_LP, [j=1:m], sum(z_w_LP[i,j] for i in 1:n) == 1)
    @constraint(model_w_LP, [i=1:n], sum(z_w_LP[i,j] for j in 1:m) <= m*x_w_LP[i])
    @objective(model_w_LP, Min, sum(f[i]*x_w_LP[i] for i in 1:n) + sum((c[i,j]*d[j]*z_w_LP[i,j] for j in 1:m) for i in 1:n))
    optimize!(model_w_LP);

    ## STRONG MODEL
    model_s = Model(Gurobi.Optimizer)
    @variable(model_s, x_s[1:n], Bin)
    @variable(model_s, z_s[1:n,1:m], Bin)
    @constraint(model_s, [j=1:m], sum(z_s[i,j] for i in 1:n) == 1)
    @constraint(model_s, [i=1:n, j=1:m], z_s[i,j] <= x_s[i])
    @objective(model_s, Min, sum(f[i]*x_s[i] for i in 1:n) + sum((c[i,j]*d[j]*z_s[i,j] for j in 1:m) for i in 1:n))
    optimize!(model_s);
    model_s_LP = Model(Gurobi.Optimizer)    ## RELAXATION
    @variable(model_s_LP, 0 <= x_s_LP[1:n] <= 1)
    @variable(model_s_LP, 0 <= z_s_LP[1:n,1:m] <= 1)
    @constraint(model_s_LP, [j=1:m], sum(z_s_LP[i,j] for i in 1:n) == 1)
    @constraint(model_s_LP, [i=1:n, j=1:m], z_s_LP[i,j] <= x_s_LP[i])
    @objective(model_s_LP, Min, sum(f[i]*x_s_LP[i] for i in 1:n) + sum((c[i,j]*d[j]*z_s_LP[i,j] for j in 1:m) for i in 1:n))
    optimize!(model_s_LP);

    return objective_value(model_w), objective_value(model_w_LP), objective_value(model_s), objective_value(model_s_LP);
end

soln = facility_loc(m, n, param_UB, param_LB, f, c, d);
# PRINT PROBLEM PARAMETERS
print("f = "*string(f)*"\n");
print("c = "*string(c)*"\n");
print("d = "*string(d)*"\n");
# PRINT OPTIMAL SOLUTIONS
print("weak_opt = "*string(soln[1])*"\n");
print("weak_opt_LP = "*string(soln[2])*"\n");
print("strong_opt = "*string(soln[3])*"\n");
print("strong_opt_LP = "*string(soln[4])*"\n");
# ABSOLUTE GAPS
print("weak_GAP_abs = "*string(soln[1]-soln[2])*"\n");
print("strong_GAP_abs = "*string(soln[3]-soln[4])*"\n");
# RELATIVE GAPS
print("weak_GAP_rel = "*string(soln[2]/soln[1])*"\n");
print("strong_GAP_rel = "*string(soln[4]/soln[3])*"\n");
