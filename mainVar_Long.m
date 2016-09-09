
clear all;
close all;
addpath(genpath('C:\Users\onur\Desktop\semivariance november\YALMIP\yalmip'))
addpath('C:/gurobi650/win64/matlab')


%set parameters
solver = 'gurobi';
alpha = .01;

delta = 0.01;                %target precision for VaR

%cols_o_ub = ceil(2*alpha*m);   %starting columns for UpperBound
eps = 1E-5;

% Result Codes
% code_lb = 0   no problem
% code_lb = 1   Infeasible or Unbounded (should not happen)
% code_lb = 2   max time reached
% code_ub = 3   Trouble, something happened

% code_in_ub = 0   no problem, feasible or infeasible or unbounded
% code_in_ub = 2   max time reached
% code_in_ub = 3   Trouble, something happened

% code_ub = 0      VaR optimality proved
% code_ub = 1      VaR optimality not proved


% code_opt = 0   no problem
% code_opt = 1   Infeasible or Unbounded (should not happen)
% code_opt = 2   max time reached
% code_opt = 3   Trouble, something happened

% Maximum Time
 Max_Time = 3600;
% Max_Time = 1;

%No_tests = 3;

No_tests = 7;

% Create results file
fileID = fopen('Test_Results.txt','w');
fclose(fileID);

     
for n=40:10:90
 for m=2000:500:3500                   
    
    fprintf('n= %g, m= %g\n', n, m);
    VaR_pos = floor(m*alpha)+1;
    load Port100_Data;
    [rowsData,colsData] = size(Port100_Data);
    Returns = Port100_Data(rowsData-m+1:rowsData,1:n);
    clear Port100_Data;

    % Compute average return info
    mean_Returns = mean(Returns);
    max_mean = max(mean_Returns);
    min_mean = min(mean_Returns);
    
    delta_mu_o = (max_mean - min_mean)/No_tests;
    
    for mu_o=min_mean+delta_mu_o:delta_mu_o:max_mean-delta_mu_o
        
        fileID = fopen('Test_Results.txt','a');
        fprintf(fileID,'%g %g %.5f ',n, m, mu_o);
        fclose(fileID);
    
        
        % Compute Big-M value
        BigM = ceil(3*(max(abs(Returns'))));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run Lower Bound Algorithm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        cols_o_lb = ceil(2*alpha*m);   %starting columns for LowerBound

        [time_lb, ind_sam_lb, var_lb, port_lb, code_lb, iter_lb]= VaR_branch_price_Long(alpha, mu_o, m,n ,Returns, mean_Returns, BigM, solver, cols_o_lb, Max_Time);
        
        fileID = fopen('Test_Results.txt','a');
        
        if code_lb == 0,
            fprintf(fileID,'%g %.5f %.5f %g ',code_lb, var_lb, time_lb, iter_lb);
        elseif code_lb >= 1
            fprintf(fileID,'%g %.5f %.5f %g ',code_lb, -101, time_lb, iter_lb);
        end;
        
        fclose(fileID);
            
        display('Lower Bound Done');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run Upper Bound Algorithm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Set VaR goal
        var_goal = var_lb+delta*abs(var_lb);

        [time_ub, code_ub, code_in_ub, iter_ub]= VaRdual_sampled_Long(alpha, var_goal ,m, n, mean_Returns, solver, Returns, mu_o, max(0, Max_Time - time_lb), delta, port_lb);

        fileID = fopen('Test_Results.txt','a');
        
        if code_ub == 0,
            fprintf(fileID,'%g %g %.5f %g ',code_ub, code_in_ub, time_ub, iter_ub);
        elseif code_ub == 1
            fprintf(fileID,'%g %g %.5f %g ',code_ub, code_in_ub, time_ub, iter_ub);
        end;
        
        fclose(fileID);
        
        display('Upper Bound Done');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Solve by evaluating MIP VaR formulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute Big-M value
        BigM = ceil(3*(max(abs(Returns'))));


        [time_opt, var_opt, port_opt, code_opt] = VaR_opt(alpha,mu_o, m,n, Returns, mean_Returns, BigM, solver, Max_Time);

        fileID = fopen('Test_Results.txt','a');
        
        if code_opt == 0,
            fprintf(fileID,'%g %.5f %.5f\n',code_opt, var_opt, time_opt);
        elseif code_opt >= 1
            fprintf(fileID,'%g %.5f %.5f\n',code_opt, -101, time_opt);
            var_opt=-101;
        end;
        
        fclose(fileID);
        
        display('IP Done');
        %table=[n, m, alpha, mu_o, var_opt, time_opt, code_opt, 
        table=[n, m, alpha, mu_o, code_opt, var_opt, time_opt,code_lb, var_lb, time_lb, iter_lb, code_ub, code_in_ub, time_ub, iter_ub];
        dlmwrite('6-14-VaRResultsyalmipalpha01.csv',table,'-append')

    end;
 end;
end;



