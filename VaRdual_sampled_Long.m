function [time, code, code_in, iter]= VaRdual_sampled_Long(alpha, VaR, m, n, mean_Returns, solver, Returns, mu_o, Max_Time, delta, port_lb);


%% start with no success

code = 1;

m_orig = size(Returns,1);

m1 = floor(alpha*m);

beta = 0.01;

added_samples = ceil(beta*m1)+1;

%% start clock & iterations

tic;
iter = 1;
code_in = 0;

Port_returns_lb = Returns*port_lb;
[Port_returns_lb_s, ind_lb_s] = sort(Port_returns_lb);

ind_ub_s = ind_lb_s;

I = [ind_ub_s(1:floor((1+2*beta)*m1))];

while code && code_in == 0
    
    %Returns_ub = Returns(ind_ub_s(1:floor((1+beta*iter)*m1)+2*iter),:);
    
    Returns_ub = Returns(I,:);
    
    BigM = ceil(3*(max(abs(Returns_ub'))));
    
    m = size(Returns_ub,1);

    % Variables (sample_diferences(y), portfolio(x), VaR(v))
    y = binvar(m,1);
    x = sdpvar(n,1);
    
    % Objective
    obj = mean_Returns*x;
    
    % Constraints
    F = [];
    
    % Sum of positions = 1
    F = F+[sum(x) == 1];
    
    % VaR  constraints
    F = F + [sum(y) <= m1];
    
    % Big M to the left constraint
    F = F + [ BigM'.*y >=  VaR - Returns_ub*x];
    
    % Upper and Lower bounds
    F = F + [0 <= x <= 1];
    
    % Solve Problem
    diagnostics = solvesdp(F,-obj,sdpsettings('solver',solver,'gurobi.TimeLimit',Max_Time,'verbose',0));
    
    Port_returns_ub = Returns*double(x);
    [Port_returns_ub_s, ind_ub_s] = sort(Port_returns_ub);
    
    %[ind_lb_s'; Port_returns_lb_s'; ind_ub_s'; Port_returns_ub_s']

    %diagnostics.problem
    

    if diagnostics.problem == 0
        code_in = 0;
    elseif diagnostics.problem == 1 || diagnostics.problem == 12  %It cannot be Unbounded!? 12: Infeasible or Unbounded
         code_in = 0;
    elseif diagnostics.problem == 3
        code_in = 2;
    else
        display(yalmiperror(diagnostics.problem))
        code_in = 3;
    end;
    
    time_ub = toc;
    
    if double(obj) <= mu_o || diagnostics.problem == 1 || diagnostics.problem == 12  % second term is "infeasibility"
        code = 0;
    end;
    
    if time_ub >= Max_Time
        code=0;
        code_in=2;
    end
    
%     % two of lower bound
    J = setdiff(ind_lb_s,I);
    Port_returns_J = Returns(J,:)*port_lb;
    [Port_returns_J_s, ind_J_s] = sort(Port_returns_J);
    
    I = [I;J(ind_J_s(1:added_samples))];
%     
%     % two of upper bound
%     J = setdiff(ind_lb_s,I);
%     Port_returns_J = Returns(J,:)*double(x);
%     [Port_returns_J_s, ind_J_s] = sort(Port_returns_J);
%     
%     I = [I;J(ind_J_s(1:4))];
    
    % one from each
%     J = setdiff(ind_lb_s,I);
%     Port_returns_J = Returns(J,:)*port_lb;
%     [Port_returns_J_s, ind_J_s] = sort(Port_returns_J);
%     
%      I = [I;J(ind_J_s(1:2))];
%     
%     J = setdiff(ind_lb_s,I);
%     Port_returns_J = Returns(J,:)*double(x);
%     [Port_returns_J_s, ind_J_s] = sort(Port_returns_J);
    
  %  I = [I;J(ind_J_s(1:2))];
    
    
    iter = iter+1;
end;

time = toc;


