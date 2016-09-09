function [time, var, port, code]= VaR_opt(alpha, mu_o, m,n ,Returns, mean_Returns, BigM, solver, Max_Time);


tic;

% Variables (sample_indicators(y), VaR(v), portfolio(x))
y = binvar(m,1);
v = sdpvar(1,1);
x = sdpvar(n,1);

% Objective
obj = v;

% Constraints
F = [];

% Sum of positions = 1
F = F+[sum(x) == 1];

% Minimum expected return constraint

F = F + [mean_Returns*x >= mu_o];
% Sum of y's constraint

F = F + [sum(y) == floor(m*alpha)];

% Big M to the left constraint
F = F + [ BigM'.*y >=  v - Returns*x];

% Upper and Lower bounds
F = F + [min(min(Returns))<= v <= max(max(Returns))];
F = F + [0 <= x <= 1];

% Solve Problem
diagnostics = solvesdp(F,-obj,sdpsettings('solver',solver,'gurobi.TimeLimit',Max_Time,'verbose',0));


 if diagnostics.problem == 0
        code = 0;
    elseif diagnostics.problem == 1 || diagnostics.problem == 12  %It cannot be Unbounded!? 12: Infeasible or Unbounded
         code = 1;
    elseif diagnostics.problem == 3
        code = 2;
    else
        code = 3;
    end;


% set output
time = toc;
port = double(x);
ind_sam = double(y);
var = double(v);
