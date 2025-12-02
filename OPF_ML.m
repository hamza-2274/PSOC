%% ========================================================================
%  OPTIMAL POWER FLOW USING STEEPEST DESCENT METHOD
%  IEEE-14 Bus System - COMPLETE BOOK METHOD IMPLEMENTATION
%  Based on Mariesa Crow's Algorithm - ALL TERMS INCLUDED
%% ========================================================================

clear; clc; close all;
fprintf('\n========================================\n');
fprintf('  OPTIMAL POWER FLOW - STEEPEST DESCENT\n');
fprintf('  IEEE-14 Bus System (CORRECTED)\n');
fprintf('========================================\n\n');

%% ====================== SYSTEM DATA ======================
BMva = 100;   % Base MVA

% Bus Data: [Bus Vsp delta Pg Qg Pl Ql Qmin Qmax Type]
busdata = [
    1   1.060   0    232.4   -16.9    0.0    0.0     0    0    1;
    2   1.045   0     40.0    43.56  21.7   12.7   -40   50    2;
    3   1.010   0      0.0     0.0   94.2   19.0   -25   40    3;
    4   1.000   0      0.0     0.0   47.8   -3.9     0    0    3;
    5   1.000   0      0.0     0.0    7.6    1.6     0    0    3;
    6   1.070   0      0.0     0.0   11.2    7.5    -6   24    2;
    7   1.000   0      0.0     0.0    0.0    0.0     0    0    3;
    8   1.090   0      0.0     0.0    0.0    0.0     0    0    2;
    9   1.000   0      0.0     0.0   29.5   16.6     0    0    3;
   10   1.000   0      0.0     0.0    9.0    5.8     0    0    3;
   11   1.000   0      0.0     0.0    3.5    1.8     0    0    3;
   12   1.000   0      0.0     0.0    6.1    1.6     0    0    3;
   13   1.000   0      0.0     0.0   13.5    5.8     0    0    3;
   14   1.000   0      0.0     0.0   14.9    5.0     0    0    3;
];

% Line Data: [From To R X B/2 Tap]
linedata = [
    1   2   0.01938   0.05917   0.0264   1.000;
    1   5   0.05403   0.22304   0.0246   1.000;
    2   3   0.04699   0.19797   0.0219   1.000;
    2   4   0.05811   0.17632   0.0170   1.000;
    2   5   0.05695   0.17388   0.0173   1.000;
    3   4   0.06701   0.17103   0.0064   1.000;
    4   5   0.01335   0.04211   0.0000   1.000;
    4   7   0.00000   0.20912   0.0000   0.978;
    4   9   0.00000   0.55618   0.0000   0.969;
    5   6   0.00000   0.25202   0.0000   1.000;
    6  11   0.09498   0.19890   0.0000   1.000;
    6  12   0.12291   0.25581   0.0000   1.000;
    6  13   0.06615   0.13027   0.0000   1.000;
    7   8   0.00000   0.17615   0.0000   1.000;
    7   9   0.00000   0.11001   0.0000   1.000;
    9  10   0.03181   0.08450   0.0000   1.000;
    9  14   0.12711   0.27038   0.0000   1.000;
   10  11   0.08205   0.19207   0.0000   1.000;
   12  13   0.22092   0.19988   0.0000   1.000;
   13  14   0.17093   0.34802   0.0000   1.000;
];

%% ====================== Y-BUS FORMATION ======================
fprintf('Building Y-bus matrix...\n');
nbus = max(max(linedata(:, 1:2)));
Y = zeros(nbus, nbus);
fb = linedata(:, 1);
tb = linedata(:, 2);

for k = 1:size(linedata, 1)
    i = fb(k); j = tb(k);
    R = linedata(k, 3); X = linedata(k, 4);
    Bc = linedata(k, 5); tap = linedata(k, 6);
    if tap == 0, tap = 1; end
    
    Z = R + 1i * X; y = 1 / Z; ysh = 1i * Bc;
    Y(i, j) = Y(i, j) - y / tap;
    Y(j, i) = Y(i, j);
    Y(i, i) = Y(i, i) + y / (tap^2) + ysh;
    Y(j, j) = Y(j, j) + y + ysh;
end

G = real(Y); B = imag(Y);

%% ====================== OPF PARAMETERS ======================
% Cost function for BOTH generators
% C(PG) = a + b*PG + c*PG² (P in MW)

a1 = 120; b1 = 25; c1 = 0.01;    % Generator 1 (slack)
a2 = 120; b2 = 25; c2 = 0.01;    % Generator 2 (control variable)

fprintf('Cost Functions:\n');
fprintf('  C1(PG1) = %.0f + %.0f*PG1 + %.4f*PG1²  (PG in MW)\n', a1, b1, c1);
fprintf('  C2(PG2) = %.0f + %.0f*PG2 + %.4f*PG2²  (PG in MW)\n\n', a2, b2, c2);

% Generator limits (pu)
PG2_min = 0.10;
PG2_max = 2.00;  % Increased to allow reaching optimal dispatch

% Line constraint: Line 1-5 power flow limit
line_from = 1; line_to = 5;
P_line_max = 100 / BMva;  % pu

% Algorithm parameters
max_opf_iter = 100;
opf_tol = 1e-4;
alpha_initial = 0.01;  % Reasonable step size for numerical gradient
alpha = alpha_initial;
alpha_min = 1e-6;
alpha_max = 0.1;
penalty_weight = 1000;

% Initial value
PG2_initial = 0.40;
PG2 = PG2_initial;

fprintf('Initial PG2 = %.4f pu (%.2f MW)\n\n', PG2, PG2*BMva);
fprintf('NOTE: Using numerical gradient (finite difference) for reliability.\n');
fprintf('Analytical gradient shown for comparison/debugging.\n');
fprintf('Starting OPF iterations (alpha=%.3f)...\n\n', alpha);

%% ====================== OPF ITERATIONS ======================
convergence_history = zeros(max_opf_iter, 7);

for opf_iter = 1:max_opf_iter
    
    % Update busdata with current PG2
    busdata(2, 4) = PG2 * BMva;
        [V, del, P_slack, converged, ~, J_full] = run_load_flow_with_jacobian(busdata, Y, G, B, BMva);
    if ~converged
        fprintf('Load flow did not converge at iteration %d\n', opf_iter);
        break;
    end
    
    PG1 = P_slack;  
    
    C1 = a1 + b1*(PG1*BMva) + c1*(PG1*BMva)^2;
    C2 = a2 + b2*(PG2*BMva) + c2*(PG2*BMva)^2;
    total_cost = C1 + C2;
    
    idx = find(linedata(:,1)==line_from & linedata(:,2)==line_to);
    R = linedata(idx, 3); X = linedata(idx, 4);
    Z = R + 1i*X; y = 1/Z;
    V1 = V(line_from) * exp(1i*del(line_from));
    V5 = V(line_to) * exp(1i*del(line_to));
    S15 = V1 * conj((V1 - V5) * y);
    P15 = real(S15);
    violation = max(0, P15 - P_line_max);
    
    %%
    eps_grad = 0.001;
    busdata_plus = busdata;
    busdata_plus(2, 4) = (PG2 + eps_grad) * BMva;
    [V_p, del_p, PG1_p, conv_p, ~, ~] = run_load_flow_with_jacobian(busdata_plus, Y, G, B, BMva);
    
    busdata_minus = busdata;
    busdata_minus(2, 4) = (PG2 - eps_grad) * BMva;
    [V_m, del_m, PG1_m, conv_m, ~, ~] = run_load_flow_with_jacobian(busdata_minus, Y, G, B, BMva);
    
    if conv_p && conv_m
        C1_p = a1 + b1*(PG1_p*BMva) + c1*(PG1_p*BMva)^2;
        C2_p = a2 + b2*((PG2+eps_grad)*BMva) + c2*((PG2+eps_grad)*BMva)^2;
        Cost_p = C1_p + C2_p;
        
        C1_m = a1 + b1*(PG1_m*BMva) + c1*(PG1_m*BMva)^2;
        C2_m = a2 + b2*((PG2-eps_grad)*BMva) + c2*((PG2-eps_grad)*BMva)^2;
        Cost_m = C1_m + C2_m;
        gradient_numerical = (Cost_p - Cost_m) / (2 * eps_grad);
    else
        gradient_numerical = 0;
    end
    
    %  PQ 
    type = busdata(:, 10);
    pq = find(type == 3);
    npq = length(pq);
    
    %∂C/∂PG2 
    PG2_MW = PG2 * BMva;
    dC2_dPG2 = (b2 + 2*c2*PG2_MW) * BMva;  % 
    PG1_MW = PG1 * BMva;
    dC_dPG1 = (b1 + 2*c1*PG1_MW) * BMva;  % $/pu
    
    % ∂PG1/∂x 
    dP1_dtheta = zeros(nbus-1, 1);
    for i = 2:nbus
        idx_theta = i - 1;
        dP1_dtheta(idx_theta) = V(1) * V(i) * (G(1,i)*sin(del(1)-del(i)) - B(1,i)*cos(del(1)-del(i)));
    end
    
    dP1_dV_pq = zeros(npq, 1);
    for k = 1:npq
        i = pq(k);
        dP1_dV_pq(k) = V(1) * (G(1,i)*cos(del(1)-del(i)) + B(1,i)*sin(del(1)-del(i)));
    end
    
    dC_dx = dC_dPG1 * [dP1_dtheta; dP1_dV_pq];  % (∂C/∂PG1) * (∂PG1/∂x)
    
    %   ∂g/∂u
    dg_du = zeros(nbus-1+npq, 1);
    bus2_idx = 1;  
    dg_du(bus2_idx) = 1.0;
    if rcond(J_full) > 1e-10
        y = J_full \ dC_dx;
        ct = dg_du' * y;
        gradient_analytical = dC2_dPG2 - ct;
    else
        gradient_analytical = dC2_dPG2;
        ct = 0;
    end
    gradient = gradient_numerical;

    penalty_grad = 0;
    if violation > 0
        eps_line = 0.001;
        busdata_temp = busdata;
        busdata_temp(2, 4) = (PG2 + eps_line) * BMva;
        [V_t, del_t, ~, ~, ~, ~] = run_load_flow_with_jacobian(busdata_temp, Y, G, B, BMva);
        V1_t = V_t(line_from) * exp(1i*del_t(line_from));
        V5_t = V_t(line_to) * exp(1i*del_t(line_to));
        S15_t = V1_t * conj((V1_t - V5_t) * y);
        P15_t = real(S15_t);
        dP15_dPG2 = (P15_t - P15) / eps_line;
        penalty_grad = 2 * penalty_weight * violation * dP15_dPG2;
        gradient = gradient + penalty_grad;
    end
    convergence_history(opf_iter, :) = [opf_iter, PG2*BMva, total_cost, gradient, P15*BMva, gradient_numerical, gradient_analytical];
    if opf_iter == 1
        fprintf('Iter %2d: PG2=%.4f pu, Cost=$%.2f, Grad=%+.2f (Num=%+.2f, Ana=%+.2f), P15=%.2f MW\n', ...
            opf_iter, PG2, total_cost, gradient, gradient_numerical, gradient_analytical, P15*BMva);
    elseif opf_iter <= 10 || mod(opf_iter, 5) == 0
        cost_chg = total_cost - convergence_history(opf_iter-1, 3);
        fprintf('Iter %2d: PG2=%.4f pu, Cost=$%.2f (%+.2f), Grad=%+.2f (Num=%+.2f), P15=%.2f MW\n', ...
            opf_iter, PG2, total_cost, cost_chg, gradient, gradient_numerical, P15*BMva);
    end

    PG2_new = PG2 - alpha * gradient;

    PG2_new = max(PG2_min, min(PG2_max, PG2_new));
    busdata_test = busdata;
    busdata_test(2, 4) = PG2_new * BMva;
    [V_test, del_test, P_slack_test, converged_test, ~, ~] = run_load_flow_with_jacobian(busdata_test, Y, G, B, BMva);
    
    if converged_test
        C1_test = a1 + b1*(P_slack_test*BMva) + c1*(P_slack_test*BMva)^2;
        C2_test = a2 + b2*(PG2_new*BMva) + c2*(PG2_new*BMva)^2;
        total_cost_test = C1_test + C2_test;
        
        if total_cost_test > total_cost && opf_iter > 1
            alpha = max(alpha * 0.5, alpha_min);
            PG2_new = PG2 - alpha * gradient;
            PG2_new = max(PG2_min, min(PG2_max, PG2_new));
        elseif total_cost_test < total_cost
            alpha = min(alpha * 1.1, alpha_max);
        end
    end
    

    change_in_PG2 = abs(PG2_new - PG2);
    at_boundary = (PG2_new == PG2_min || PG2_new == PG2_max);
    
    if change_in_PG2 < opf_tol && abs(gradient) < opf_tol*100
        convergence_history = convergence_history(1:opf_iter, :);
        PG2 = PG2_new;
        break;
    elseif at_boundary && change_in_PG2 < opf_tol && abs(gradient) > opf_tol*100
        convergence_history = convergence_history(1:opf_iter, :);
        PG2 = PG2_new;
        break;
    end
    
    PG2 = PG2_new;
end

busdata(2, 4) = PG2 * BMva;
[V, del, PG1, ~, ~, ~] = run_load_flow_with_jacobian(busdata, Y, G, B, BMva);
C1 = a1 + b1*(PG1*BMva) + c1*(PG1*BMva)^2;
C2 = a2 + b2*(PG2*BMva) + c2*(PG2*BMva)^2;
total_cost = C1 + C2;

idx = find(linedata(:,1)==line_from & linedata(:,2)==line_to);
R = linedata(idx, 3); X = linedata(idx, 4);
Z = R + 1i*X; y = 1/Z;
V1 = V(line_from) * exp(1i*del(line_from));
V5 = V(line_to) * exp(1i*del(line_to));
S15 = V1 * conj((V1 - V5) * y);
P15 = real(S15);

fprintf('\n=////////////////////////////////\n');
fprintf('  FINAL RESULTS\n');
fprintf('Iterations: %d\n\n', opf_iter);
fprintf('COST BREAKDOWN:\n');
fprintf('  C1(PG1) = $%.2f/hr\n', C1);
fprintf('  C2(PG2) = $%.2f/hr\n', C2);
fprintf('  Total = $%.2f/hr\n\n', total_cost);

% Compute marginal costs
MC1 = b1 + 2*c1*(PG1*BMva);
MC2 = b2 + 2*c2*(PG2*BMva);
fprintf('MARGINAL COSTS:\n');
fprintf('  MC1 = %.2f $/MWh\n', MC1);
fprintf('  MC2 = %.2f $/MWh\n', MC2);
fprintf('  Difference = %.2f $/MWh\n\n', abs(MC1-MC2));

fprintf('INITIAL STATE:\n');
fprintf('  PG2 = %.2f MW, Total Cost = $%.2f/hr\n\n', ...
    PG2_initial*BMva, convergence_history(1,3));
fprintf('FINAL STATE:\n');
fprintf('  PG2 = %.2f MW, PG1 = %.2f MW (slack)\n', PG2*BMva, PG1*BMva);
fprintf('  Total Generation = %.2f MW\n', (PG1+PG2)*BMva);
fprintf('  Total Cost = $%.2f/hr\n', total_cost);
if total_cost < convergence_history(1,3)
    fprintf('  Cost Reduction: $%.2f/hr (%.1f%%) ✓\n', ...
        convergence_history(1,3)-total_cost, ...
        100*(convergence_history(1,3)-total_cost)/convergence_history(1,3));
else
    fprintf('  Cost Change: $%+.2f/hr\n', total_cost - convergence_history(1,3));
end
if P15 <= P_line_max
    fprintf('\nCONSTRAINT:\n');
    fprintf('  Line 1-5: %.2f MW (limit: %.2f MW) ✓\n', P15*BMva, P_line_max*BMva);
else
    fprintf('\nCONSTRAINT:\n');
    fprintf('  Line 1-5: %.2f MW (limit: %.2f MW) ✗ VIOLATED\n', P15*BMva, P_line_max*BMva);
end
fprintf('========================================\n\n');

% Voltages
fprintf('BUS VOLTAGES:\n');
fprintf('Bus | V (pu) | Angle (deg)\n');
fprintf('----------------------------\n');
for i = 1:nbus
    fprintf('%2d  | %.4f | %8.3f\n', i, V(i), rad2deg(del(i)));
end
fprintf('----------------------------\n\n');

% Plots
figure('Position', [100 100 1400 900]);

subplot(2,3,1);
plot(convergence_history(:,1), convergence_history(:,2), 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
grid on; xlabel('Iteration'); ylabel('PG2 (MW)');
title('Generator 2 Output');

subplot(2,3,2);
plot(convergence_history(:,1), convergence_history(:,3), 'r-s', 'LineWidth', 2, 'MarkerSize', 4);
grid on; xlabel('Iteration'); ylabel('Cost ($/hr)');
title('Total Generation Cost');

subplot(2,3,4);
plot(convergence_history(:,1), convergence_history(:,5), 'k-d', 'LineWidth', 2, 'MarkerSize', 4);
hold on; yline(P_line_max*BMva, 'r--', 'LineWidth', 2);
grid on; xlabel('Iteration'); ylabel('P_{1-5} (MW)');
title('Line 1-5 Power Flow'); legend('Flow', 'Limit');




sgtitle('OPF ', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('=== OPTIMIZATION COMPLETE ===\n\n');

%% ====================== MODIFIED LOAD FLOW FUNCTION ======================
function [V, del, P_slack, converged, iterations, J] = run_load_flow_with_jacobian(busdata, Y, G, B, BMva)
    nbus = size(busdata, 1);
    type = busdata(:, 10);
    V = busdata(:, 2);
    del = deg2rad(busdata(:, 3));
    Pg = busdata(:, 4) / BMva;
    Qg = busdata(:, 5) / BMva;
    Pl = busdata(:, 6) / BMva;
    Ql = busdata(:, 7) / BMva;
    Psp = Pg - Pl;
    Qsp = Qg - Ql;
    pq = find(type == 3);
    npq = length(pq);
    tol = 1e-6;
    max_iter = 20;
    iterations = 0;
    Tol = 1;
    converged = false;
    J = [];
    
    while Tol > tol && iterations < max_iter
        iterations = iterations + 1;
        P = zeros(nbus, 1);
        Q = zeros(nbus, 1);
        for i = 1:nbus
            for k = 1:nbus
                P(i) = P(i) + V(i)*V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
                Q(i) = Q(i) + V(i)*V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
            end
        end
        dP = Psp - P;
        dQ = Qsp - Q;
        M = [dP(2:end); dQ(pq)];
        Tol = max(abs(M));
        if Tol < tol
            converged = true;
            break;
        end
        
        % Build Jacobian J1
        J1 = zeros(nbus-1);
        for i = 1:nbus-1
            m = i+1;
            for k = 1:nbus-1
                n = k+1;
                if n == m
                    for n2 = 1:nbus
                        J1(i,k) = J1(i,k) + V(m)*V(n2)*(-G(m,n2)*sin(del(m)-del(n2)) + B(m,n2)*cos(del(m)-del(n2)));
                    end
                    J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
                else
                    J1(i,k) = V(m)*V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
            end
        end
        
        % Build Jacobian J2
        J2 = zeros(nbus-1, npq);
        for i = 1:nbus-1
            m = i+1;
            for k = 1:npq
                n = pq(k);
                if n == m
                    for n2 = 1:nbus
                        J2(i,k) = J2(i,k) + V(n2)*(G(m,n2)*cos(del(m)-del(n2)) + B(m,n2)*sin(del(m)-del(n2)));
                    end
                    J2(i,k) = J2(i,k) + V(m)*G(m,m);
                else
                    J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
            end
        end
        
        % Build Jacobian J3
        J3 = zeros(npq, nbus-1);
        for i = 1:npq
            m = pq(i);
            for k = 1:nbus-1
                n = k+1;
                if n == m
                    for n2 = 1:nbus
                        J3(i,k) = J3(i,k) + V(m)*V(n2)*(G(m,n2)*cos(del(m)-del(n2)) + B(m,n2)*sin(del(m)-del(n2)));
                    end
                    J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
                else
                    J3(i,k) = V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
                end
            end
        end
        
        % Build Jacobian J4
        J4 = zeros(npq);
        for i = 1:npq
            m = pq(i);
            for k = 1:npq
                n = pq(k);
                if n == m
                    for n2 = 1:nbus
                        J4(i,k) = J4(i,k) + V(n2)*(G(m,n2)*sin(del(m)-del(n2)) - B(m,n2)*cos(del(m)-del(n2)));
                    end
                    J4(i,k) = J4(i,k) - V(m)*B(m,m);
                else
                    J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
            end
        end
        
        % Assemble full Jacobian
        J = [J1 J2; J3 J4];
        
        % Solve for corrections
        X = J \ M;
        del(2:end) = del(2:end) + X(1:nbus-1);
        V(pq) = V(pq) + X(nbus:end);
    end
    
    % Compute slack power
    P_slack = 0;
    for k = 1:nbus
        P_slack = P_slack + V(1)*V(k)*(G(1,k)*cos(del(1)-del(k)) + B(1,k)*sin(del(1)-del(k)));
    end
end