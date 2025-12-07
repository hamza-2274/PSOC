%% ========================================================================
%  ECONOMIC DISPATCH - LOSSLESS (λ-ITERATION WITH LIMITS)

clear; clc; close all;
PLoad = [0; 0.217; 0.942; 0.478; 0.076; 0.112; 0; 0; 0.295; 0.09; 0.035; 0.061; 0.135; 0.149];
PD = sum(PLoad);  
fprintf('Total System Demand:\n');
fprintf('  PD = %.4f pu (%.2f MW at 100 MVA base)\n\n', PD, PD*100);
a_MW = [100; 120; 150; 180; 200];          
b_MW = [20; 25; 30; 35; 40];                
c_MW = [0.005; 0.010; 0.012; 0.015; 0.020]; 
a = a_MW; b = b_MW * 100;     c = c_MW * 10000;    
Pmin_MW = [50; 20; 15; 10; 10];
Pmax_MW = [200; 150; 180; 200; 200];
Pmin = Pmin_MW / 100;Pmax = Pmax_MW / 100;
n = length(a);  

% Initialize
lambda = 2500;    
max_iter = 100;
tol = 1e-6;
alpha_step = 50;   


fprintf('Iter | Lambda ($/MWh) | ΣPG (pu) | ΣPG (MW) | Error (pu)\n');
fprintf('-----+----------------+----------+----------+-----------\n');

for iter = 1:max_iter
    % C/P = lamda: b + 2cP = λ ->  P = 
    PG = (lambda - b) ./ (2*c);
    PG = max(Pmin, min(Pmax, PG));
    PG_total = sum(PG);
    error = PG_total - PD;
    
    if iter <= 5 || mod(iter, 10) == 0 || abs(error) < tol
        fprintf('%4d | %13.2f | %8.4f | %8.2f | %+9.6f\n', ...
            iter, lambda, PG_total, PG_total*100, error);
    end
    
    % Check convergence
    if abs(error) < tol
        fprintf('\n conv at itera %d\n', iter);
        break;
    end

    lambda = lambda - alpha_step * error;
end


%% ====================== CALCULATE FINAL RESULTS ======================
% Calculate costs and incremental costs
C = zeros(n,1);
IC = zeros(n,1);

for i = 1:n
    C(i) = a(i) + b(i)*PG(i) + c(i)*PG(i)^2;
    IC(i) = b(i) + 2*c(i)*PG(i);
end

C_total = sum(C);

%% ====================== DISPLAY RESULTS ======================
fprintf('========================================\n');
fprintf('  ECONOMIC DISPATCH RESULTS\n');
fprintf('========================================\n\n');

fprintf('GENERATOR DISPATCH:\n');
fprintf('--------------------------------------------------------------------------------\n');
fprintf('Gen | Output (pu) | Output (MW) | IC ($/MWh) | Cost ($/hr) | Status\n');
fprintf('--------------------------------------------------------------------------------\n');
for i = 1:n
    % Check if at limits
    status = '';
    if abs(PG(i) - Pmin(i)) < 0.001
        status = '[at Pmin]';
    elseif abs(PG(i) - Pmax(i)) < 0.001
        status = '[at Pmax]';
    end
    
    fprintf(' G%d | %10.4f | %10.2f | %10.2f | %10.2f | %s\n', ...
        i, PG(i), PG(i)*100, IC(i), C(i), status);
end
fprintf('--------------------------------------------------------------------------------\n');
fprintf('TOTAL| %10.4f | %10.2f |            | %10.2f |\n', ...
    sum(PG), sum(PG)*100, C_total);
fprintf('--------------------------------------------------------------------------------\n\n');

fprintf('SUMMARY:\n');
fprintf('  System Lambda (λ) = %.2f $/MWh\n', lambda);
fprintf('  Total Generation = %.4f pu (%.2f MW)\n', sum(PG), sum(PG)*100);
fprintf('  Total Load = %.4f pu (%.2f MW)\n', PD, PD*100);
fprintf('  Power Balance = %.6f pu (%.4f MW)\n', sum(PG)-PD, (sum(PG)-PD)*100);
fprintf('  Total Cost = $%.2f/hr\n\n', C_total);
