function [time, ind_sam, var, port, code, iter]= VaR_branch_price_Long(alpha, mu_o, m,n ,Returns, mean_Returns, BigM, solver, cols_o_lb, Max_Time);
% cols_o_lb: Starting columns for Lower Bound VaR algorithm

tic;

% Start samples to be allowed to be binary
Iy = (1:cols_o_lb)';

non_STOP = 1;
iter = 0;

time_lb = 0;
code_in = 0;

while time_lb < Max_Time  && (non_STOP) && code_in == 0

    % Variables (sample_indicators(y), VaR(v), portfolio(x))
    y = binvar(m,1);
    v = sdpvar(1,1);
    x = sdpvar(n,1);
    
    y(setdiff(1:m,Iy)) = 0;
   
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
    
    % Set y's to zero
    %F = F + [y(setdiff(1:m,Iy)) == 0];
    
    % Big M to the left constraint
    F = F + [ BigM'.*y >=  v - Returns*x];
    
    % Upper and Lower bounds
    F = F + [min(min(Returns))<= v <= max(max(Returns))];
    F = F + [0 <= x <= 1];

    % Solve Problem
    diagnostics = solvesdp(F,-obj,sdpsettings('solver',solver,'gurobi.TimeLimit',Max_Time,'verbose',0));
    
    if diagnostics.problem == 0
        code_in = 0;
    elseif diagnostics.problem == 1 || diagnostics.problem == 12  %It cannot be Unbounded!? 12: Infeasible or Unbounded
          code_in = 1;
    elseif diagnostics.problem == 3
        code_in = 2;
    else
        code_in = 3;
    end
    
    %% gurobi.TimeLimit
    
    % Update samples to be allowed to be binary
    Istar =  find(double(y) > 0.99);

    slack = double(v - Returns*x);
    
    Islack = find(slack >-1e-3);

    oldIy = Iy;
    
    Iy = [Istar; Islack];

    Iy = sort(Iy);

    if length(oldIy) == length(Iy)
        non_STOP = sum(oldIy ~= Iy);
    end
    %Istar'
    %Inew
    
    iter = iter+1;
    time_lb=toc;
    if time_lb >= Max_Time
        code_in=2;
    end
end

% set output
time = toc
port = double(x);
ind_sam = double(y);
var = double(v)

code= code_in;

