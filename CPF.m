%% ========================================================================
%  CONTINUOUS POWER FLOW (CPF) ANALYSIS - IEEE-14 Bus System
%  
%  This program implements Continuation Power Flow method to trace PV
%  curves and analyze voltage stability in the stable region.
%  ========================================================================

clear; clc; close all;
fprintf('\n==== CONTINUOUS POWER FLOW (CPF) - IEEE-14 BUS ====\n\n');

%% ====================== SYSTEM DATA ======================
BMva = 100;   % Base MVA

% Bus Data Format:
% Columns: Bus | Vsp | delta | Pg | Qg | Pl | Ql | Qmin | Qmax | Type
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

% Line Data
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
    i = fb(k);
    j = tb(k);
    R = linedata(k, 3);
    X = linedata(k, 4);
    Bc = linedata(k, 5);
    tap = linedata(k, 6);
    if tap == 0, tap = 1; end
    
    Z = R + 1i * X;
    y = 1 / Z;
    ysh = 1i * Bc;
    
    Y(i, j) = Y(i, j) - y / tap;
    Y(j, i) = Y(i, j);
    Y(i, i) = Y(i, i) + y / (tap^2) + ysh;
    Y(j, j) = Y(j, j) + y + ysh;
end

G = real(Y);
B = imag(Y);

%% ====================== BASE CASE PARAMETERS ======================
type = busdata(:, 10);
V0 = busdata(:, 2);
del0 = deg2rad(busdata(:, 3));
Pg0 = busdata(:, 4) / BMva;
Qg0 = busdata(:, 5) / BMva;
Pl0 = busdata(:, 6) / BMva;  % Base load P
Ql0 = busdata(:, 7) / BMva;  % Base load Q

pq = find(type == 3);
pv = find(type == 2);
npq = length(pq);
npv = length(pv);

%% ====================== CPF INITIALIZATION ======================
fprintf('Initializing CPF at lambda = 0...\n\n');

% CPF Parameters
lambda = 0;           % Loading parameter
lambda_max = 5.0;     % Maximum lambda to explore
d_lambda = 0.05;      % Initial step size
tol = 1e-6;           % Convergence tolerance
max_iter = 20;        % Max NR iterations per step

% Storage for PV curves
lambda_trace = [];
V_trace = [];

% Current state
V = V0;
del = del0;

% Direction of load increase (uniform across all loads)
Pl_dir = Pl0;  % Direction vector for P
Ql_dir = Ql0;  % Direction vector for Q

fprintf('====================================================\n');
fprintf('  CPF SETTINGS\n');
fprintf('====================================================\n');
fprintf('  Lambda range: [0, %.2f]\n', lambda_max);
fprintf('  Initial step: %.4f\n', d_lambda);
fprintf('  Convergence tolerance: %.2e\n', tol);
fprintf('====================================================\n\n');

%% ====================== STEP 1: SOLVE BASE CASE (λ=0) ======================
fprintf('Step 1: Solving base case at lambda = 0...\n');

Pl = Pl0;
Ql = Ql0;
Psp = Pg0 - Pl;
Qsp = Qg0 - Ql;

[V, del, converged, iter] = solve_power_flow(V, del, Psp, Qsp, Y, G, B, ...
                                              type, pq, nbus, tol, max_iter);

if ~converged
    error('Base case power flow did not converge!');
end

fprintf('  Base case converged in %d iterations\n', iter);
fprintf('  Storing initial point (lambda = %.4f)\n\n', lambda);

lambda_trace = [lambda_trace; lambda];
V_trace = [V_trace, V];

%% ====================== CPF CONTINUATION LOOP ======================
fprintf('Starting CPF continuation...\n');
fprintf('----------------------------------------------\n');

step_count = 0;
cpf_failed = false;

while lambda < lambda_max && ~cpf_failed
    step_count = step_count + 1;
    
    %% PREDICTOR STEP
    % Compute tangent vector using augmented Jacobian
    
    % Update load specifications
    Pl = Pl0 + lambda * Pl_dir;
    Ql = Ql0 + lambda * Ql_dir;
    Psp = Pg0 - Pl;
    Qsp = Qg0 - Ql;
    
    % Compute current Jacobian
    J = compute_jacobian(V, del, G, B, nbus, pq);
    
    % Augment Jacobian with load direction
    % dP/dlambda = -Pl_dir, dQ/dlambda = -Ql_dir
    dP_dlambda = -Pl_dir(2:end);
    dQ_dlambda = -Ql_dir(pq);
    
    % Augmented system: [J | [dP_dlambda; dQ_dlambda]] * [dx; dlambda]' = 0
    % with normalization constraint
    
    % Tangent vector calculation
    F_lambda = [dP_dlambda; dQ_dlambda];
    
    % Solve: J * dx = F_lambda to get tangent direction
    dx = J \ F_lambda;
    
    % Normalize tangent vector
    tangent = [dx; 1];  % [dTheta; dV; dlambda]
    tangent = tangent / norm(tangent);
    
    % Predictor: estimate next point
    d_theta_pred = tangent(1:nbus-1) * d_lambda;
    d_V_pred = tangent(nbus:end-1) * d_lambda;
    d_lambda_pred = tangent(end) * d_lambda;
    
    % Predicted state
    del_pred = del;
    del_pred(2:end) = del(2:end) + d_theta_pred;
    V_pred = V;
    V_pred(pq) = V(pq) + d_V_pred;
    lambda_pred = lambda + d_lambda_pred;
    
    %% CORRECTOR STEP
    % Use Newton-Raphson to correct the predicted state
    
    Pl = Pl0 + lambda_pred * Pl_dir;
    Ql = Ql0 + lambda_pred * Ql_dir;
    Psp = Pg0 - Pl;
    Qsp = Qg0 - Ql;
    
    [V_corr, del_corr, converged, iter] = solve_power_flow(V_pred, del_pred, ...
                                          Psp, Qsp, Y, G, B, type, pq, nbus, tol, max_iter);
    
    if ~converged
        fprintf('  Step %d: Corrector failed at lambda = %.4f\n', step_count, lambda_pred);
        fprintf('  Reducing step size and retrying...\n');
        d_lambda = d_lambda * 0.5;
        
        if d_lambda < 1e-4
            fprintf('  Step size too small. Stopping CPF.\n');
            cpf_failed = true;
        end
        continue;
    end
    
    % Update state
    V = V_corr;
    del = del_corr;
    lambda = lambda_pred;
    
    % Store results
    lambda_trace = [lambda_trace; lambda];
    V_trace = [V_trace, V];
    
    % Display progress
    if mod(step_count, 5) == 0
        fprintf('  Step %3d: lambda = %.4f, Min V = %.4f pu (Bus %d)\n', ...
                step_count, lambda, min(V), find(V == min(V), 1));
    end
    
    % Adaptive step size control
    if iter <= 3
        d_lambda = min(d_lambda * 1.2, 0.1);  % Increase step
    elseif iter >= 6
        d_lambda = d_lambda * 0.8;  % Decrease step
    end
end

fprintf('----------------------------------------------\n');
fprintf('CPF completed: %d continuation steps\n', step_count);
fprintf('Final lambda: %.4f\n\n', lambda);

%% ====================== IDENTIFY DIVERGENCE POINT ======================
% Find the critical loading point (where stable region ends)
lambda_critical = lambda;
V_critical = V;

% Find voltage at critical buses
bus_14_critical_V = V_trace(14, end);
bus_3_critical_V = V_trace(3, end);

fprintf('====================================================\n');
fprintf('  CPF ANALYSIS COMPLETE\n');
fprintf('====================================================\n');
fprintf('  Total continuation steps: %d\n', step_count);
fprintf('  Lambda range covered: [0, %.4f]\n', lambda);
fprintf('  Stable region explored successfully\n');
fprintf('====================================================\n\n');

fprintf('====================================================\n');
fprintf('  CRITICAL POINT IDENTIFICATION\n');
fprintf('  (End of Stable Region - Part 1 of PV Curve)\n');
fprintf('====================================================\n');
fprintf('  Critical Loading Parameter (λ_critical): %.4f\n', lambda_critical);
fprintf('  Load Increase Factor: %.2f%% of base load\n', (1 + lambda_critical) * 100);
fprintf('  \n');
fprintf('  VOLTAGES AT CRITICAL POINT:\n');
fprintf('  --------------------------------------------------\n');
fprintf('    Bus 14 (Load Bus): %.4f pu\n', bus_14_critical_V);
fprintf('    Bus 3  (Weakest):  %.4f pu\n', bus_3_critical_V);
fprintf('    System Minimum:    %.4f pu at Bus %d\n', min(V_critical), find(V_critical == min(V_critical), 1));
fprintf('  --------------------------------------------------\n');
fprintf('  \n');
fprintf('  INTERPRETATION:\n');
fprintf('  - Beyond λ = %.4f, the power flow diverges\n', lambda_critical);
fprintf('  - This marks the "nose point" of the PV curve\n');
fprintf('  - System can support up to %.2f%% load increase\n', lambda_critical * 100);
fprintf('  - Voltage stability margin exhausted\n');
fprintf('====================================================\n\n');

% Critical buses (typically weakest buses in IEEE-14)
critical_buses = [4, 5, 9, 10, 14];

% Plot PV Curves
figure('Position', [100 100 1400 900]);

% Plot 1: Bus 14 PV Curve (Primary Result)
subplot(2, 3, 1);
hold on; grid on; box on;
plot(lambda_trace, V_trace(14, :), 'b-', 'LineWidth', 3);
plot(lambda_trace, V_trace(14, :), 'ro', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'r', 'MarkerIndices', 1:3:length(lambda_trace));
% Mark critical point
plot(lambda_critical, bus_14_critical_V, 'kp', 'MarkerSize', 20, ...
     'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
text(lambda_critical, bus_14_critical_V + 0.03, ...
     sprintf('Critical Point\n(λ=%.3f, V=%.3f pu)', lambda_critical, bus_14_critical_V), ...
     'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black');
xlabel('Loading Parameter λ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage Magnitude (pu)', 'FontSize', 12, 'FontWeight', 'bold');
title('PV Curve - Bus 14 (Load Bus)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0.5 1.05]);
xlim([0 lambda_critical * 1.1]);

% Plot 2: Bus 3 PV Curve (Weakest Bus)
subplot(2, 3, 2);
hold on; grid on; box on;
plot(lambda_trace, V_trace(3, :), 'r-', 'LineWidth', 3);
plot(lambda_trace, V_trace(3, :), 'bo', 'MarkerSize', 6, ...
     'MarkerFaceColor', 'b', 'MarkerIndices', 1:3:length(lambda_trace));
% Mark critical point
plot(lambda_critical, bus_3_critical_V, 'kp', 'MarkerSize', 20, ...
     'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
text(lambda_critical, bus_3_critical_V + 0.03, ...
     sprintf('Critical Point\n(λ=%.3f, V=%.3f pu)', lambda_critical, bus_3_critical_V), ...
     'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
     'BackgroundColor', 'white', 'EdgeColor', 'black');
xlabel('Loading Parameter λ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage Magnitude (pu)', 'FontSize', 12, 'FontWeight', 'bold');
title('PV Curve - Bus 3 (Weakest Bus)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0.5 1.05]);
xlim([0 lambda_critical * 1.1]);

% Plot 3: All buses
subplot(2, 3, 3);
hold on; grid on; box on;
colors = lines(nbus);
for i = 1:nbus
    plot(lambda_trace, V_trace(i, :), 'LineWidth', 1.5, 'Color', colors(i, :));
end
% Highlight Bus 14
plot(lambda_trace, V_trace(14, :), 'b-', 'LineWidth', 3);
% Mark divergence region
xline(lambda_critical, '--k', 'LineWidth', 2, 'Label', 'Divergence Point', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Loading Parameter λ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage Magnitude (pu)', 'FontSize', 12, 'FontWeight', 'bold');
title('PV Curves - All Buses', 'FontSize', 14, 'FontWeight', 'bold');
legend(arrayfun(@(x) sprintf('Bus %d', x), 1:nbus, 'UniformOutput', false), ...
       'Location', 'best', 'NumColumns', 2, 'FontSize', 8);
ylim([0.5 1.1]);

% Plot 4: Critical buses comparison
subplot(2, 3, 4);
hold on; grid on; box on;
colors_crit = lines(length(critical_buses));
for idx = 1:length(critical_buses)
    i = critical_buses(idx);
    plot(lambda_trace, V_trace(i, :), 'LineWidth', 2.5, ...
         'Color', colors_crit(idx, :), 'Marker', 'o', 'MarkerSize', 4, ...
         'MarkerIndices', 1:5:length(lambda_trace));
end
% Mark critical point
xline(lambda_critical, '--k', 'LineWidth', 2);
xlabel('Loading Parameter λ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage Magnitude (pu)', 'FontSize', 12, 'FontWeight', 'bold');
title('PV Curves - Critical Load Buses', 'FontSize', 14, 'FontWeight', 'bold');
legend(arrayfun(@(x) sprintf('Bus %d', x), critical_buses, 'UniformOutput', false), ...
       'Location', 'best', 'FontSize', 10);
ylim([0.5 1.05]);

% Plot 5: Voltage profile at different loading levels
subplot(2, 3, 5);
lambda_points = [0, 0.5, 1.0, 1.5, lambda_critical];
lambda_labels = {'λ=0 (Base)', 'λ=0.5', 'λ=1.0', 'λ=1.5', sprintf('λ=%.3f (Critical)', lambda_critical)};
hold on; grid on; box on;
for lp_idx = 1:length(lambda_points)
    lp = lambda_points(lp_idx);
    [~, idx] = min(abs(lambda_trace - lp));
    if idx <= size(V_trace, 2)
        plot(1:nbus, V_trace(:, idx), '-o', 'LineWidth', 2, ...
             'MarkerSize', 6, 'DisplayName', lambda_labels{lp_idx});
    end
end
xlabel('Bus Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage Magnitude (pu)', 'FontSize', 12, 'FontWeight', 'bold');
title('Voltage Profile at Different Loading Levels', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
xlim([1 nbus]);
ylim([0.5 1.1]);

% Plot 6: Minimum voltage vs lambda with stability margin
subplot(2, 3, 6);
V_min = min(V_trace);
hold on; grid on; box on;
plot(lambda_trace, V_min, 'r-', 'LineWidth', 2.5);
plot(lambda_trace, V_min, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', ...
     'MarkerIndices', 1:5:length(lambda_trace));
% Mark critical point
plot(lambda_critical, min(V_critical), 'kp', 'MarkerSize', 20, ...
     'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
% Shade stable region
area([0 lambda_critical], [1.1 1.1], 'FaceColor', [0.8 1 0.8], ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Stable Region');
xline(lambda_critical, '--k', 'LineWidth', 2, 'Label', ...
      sprintf('Divergence\nλ=%.3f', lambda_critical), ...
      'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Loading Parameter λ', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Minimum System Voltage (pu)', 'FontSize', 12, 'FontWeight', 'bold');
title('System Minimum Voltage vs Loading', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0.5 1.05]);
legend('Min Voltage', 'Sample Points', 'Critical Point', 'Stable Region', ...
       'Location', 'best', 'FontSize', 9);

sgtitle('Continuous Power Flow Analysis - IEEE 14-Bus System (Stable Region)', ...
        'FontSize', 16, 'FontWeight', 'bold');

%% ====================== FINAL STATISTICS ======================
fprintf('\n====================================================\n');
fprintf('  DETAILED VOLTAGE STATISTICS AT CRITICAL POINT\n');
fprintf('====================================================\n');
fprintf('At λ_critical = %.4f (End of Stable Region):\n\n', lambda_critical);

fprintf('LOAD BUS VOLTAGES:\n');
fprintf('  Bus 14: %.4f pu\n', V_trace(14, end));
fprintf('  Bus 13: %.4f pu\n', V_trace(13, end));
fprintf('  Bus 12: %.4f pu\n', V_trace(12, end));
fprintf('  Bus 11: %.4f pu\n', V_trace(11, end));
fprintf('  Bus 10: %.4f pu\n', V_trace(10, end));
fprintf('  Bus  9: %.4f pu\n', V_trace(9, end));
fprintf('  Bus  5: %.4f pu\n', V_trace(5, end));
fprintf('  Bus  4: %.4f pu\n', V_trace(4, end));
fprintf('  Bus  3: %.4f pu (WEAKEST)\n\n', V_trace(3, end));

fprintf('SYSTEM-WIDE STATISTICS:\n');
fprintf('  Minimum voltage: %.4f pu at Bus %d\n', min(V_critical), find(V_critical == min(V_critical), 1));
fprintf('  Maximum voltage: %.4f pu at Bus %d\n', max(V_critical), find(V_critical == max(V_critical), 1));
fprintf('  Average voltage: %.4f pu\n', mean(V_critical));
fprintf('  Voltage range:   %.4f pu\n\n', max(V_critical) - min(V_critical));

fprintf('LOADING INFORMATION:\n');
fprintf('  Base case total load: %.2f MW\n', sum(Pl0) * BMva);
fprintf('  Critical point load:  %.2f MW\n', sum(Pl0 + lambda_critical * Pl_dir) * BMva);
fprintf('  Load increase:        %.2f MW (%.1f%%)\n', ...
        sum(lambda_critical * Pl_dir) * BMva, lambda_critical * 100);
fprintf('  Voltage stability margin: %.1f%% load increase\n', lambda_critical * 100);
fprintf('====================================================\n\n');

fprintf('KEY FINDINGS:\n');
fprintf('----------------------------------------------\n');
fprintf('1. The PV curve for Bus 14 shows stable operation\n');
fprintf('   from λ = 0 to λ = %.4f\n', lambda_critical);
fprintf('2. Load flow diverges at λ = %.4f\n', lambda_critical);
fprintf('   (This is the "nose point" of the PV curve)\n');
fprintf('3. Bus 3 is the weakest bus with V = %.4f pu\n', bus_3_critical_V);
fprintf('4. System can handle up to %.1f%% load increase\n', lambda_critical * 100);
fprintf('5. Beyond λ_critical, Part 2 of PV curve would\n');
fprintf('   require different continuation method\n');
fprintf('----------------------------------------------\n\n');

fprintf('=== CPF SIMULATION COMPLETE ===\n\n');

%% ====================== HELPER FUNCTIONS ======================

function [V, del, converged, iter] = solve_power_flow(V, del, Psp, Qsp, Y, ...
                                     G, B, type, pq, nbus, tol, max_iter)
    % Newton-Raphson power flow solver
    converged = false;
    iter = 0;
    
    while iter < max_iter
        iter = iter + 1;
        
        % Calculate power injections
        P = zeros(nbus, 1);
        Q = zeros(nbus, 1);
        
        for i = 1:nbus
            for k = 1:nbus
                P(i) = P(i) + V(i) * V(k) * (G(i,k) * cos(del(i) - del(k)) + ...
                                              B(i,k) * sin(del(i) - del(k)));
                Q(i) = Q(i) + V(i) * V(k) * (G(i,k) * sin(del(i) - del(k)) - ...
                                              B(i,k) * cos(del(i) - del(k)));
            end
        end
        
        % Calculate mismatches
        dP = Psp - P;
        dQ = Qsp - Q;
        M = [dP(2:end); dQ(pq)];
        Tol = max(abs(M));
        
        if Tol < tol
            converged = true;
            break;
        end
        
        % Compute Jacobian
        J = compute_jacobian(V, del, G, B, nbus, pq);
        
        % Solve for corrections
        X = J \ M;
        
        % Update state variables
        del(2:end) = del(2:end) + X(1:nbus-1);
        V(pq) = V(pq) + X(nbus:end);
    end
end

function J = compute_jacobian(V, del, G, B, nbus, pq)
    % Compute Jacobian matrix for power flow
    npq = length(pq);
    
    % J1: dP/d(delta)
    J1 = zeros(nbus - 1);
    for i = 1:nbus - 1
        m = i + 1;
        for k = 1:nbus - 1
            n = k + 1;
            if n == m
                for n2 = 1:nbus
                    J1(i, k) = J1(i, k) + V(m) * V(n2) * ...
                              (-G(m,n2) * sin(del(m) - del(n2)) + ...
                               B(m,n2) * cos(del(m) - del(n2)));
                end
                J1(i, k) = J1(i, k) - V(m)^2 * B(m, m);
            else
                J1(i, k) = V(m) * V(n) * ...
                          (G(m,n) * sin(del(m) - del(n)) - ...
                           B(m,n) * cos(del(m) - del(n)));
            end
        end
    end
    
    % J2: dP/dV
    J2 = zeros(nbus - 1, npq);
    for i = 1:nbus - 1
        m = i + 1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n2 = 1:nbus
                    J2(i, k) = J2(i, k) + V(n2) * ...
                              (G(m,n2) * cos(del(m) - del(n2)) + ...
                               B(m,n2) * sin(del(m) - del(n2)));
                end
                J2(i, k) = J2(i, k) + V(m) * G(m, m);
            else
                J2(i, k) = V(m) * ...
                          (G(m,n) * cos(del(m) - del(n)) + ...
                           B(m,n) * sin(del(m) - del(n)));
            end
        end
    end
    
    % J3: dQ/d(delta)
    J3 = zeros(npq, nbus - 1);
    for i = 1:npq
        m = pq(i);
        for k = 1:nbus - 1
            n = k + 1;
            if n == m
                for n2 = 1:nbus
                    J3(i, k) = J3(i, k) + V(m) * V(n2) * ...
                              (G(m,n2) * cos(del(m) - del(n2)) + ...
                               B(m,n2) * sin(del(m) - del(n2)));
                end
                J3(i, k) = J3(i, k) - V(m)^2 * G(m, m);
            else
                J3(i, k) = V(m) * V(n) * ...
                          (-G(m,n) * cos(del(m) - del(n)) - ...
                           B(m,n) * sin(del(m) - del(n)));
            end
        end
    end
    
    % J4: dQ/dV
    J4 = zeros(npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n2 = 1:nbus
                    J4(i, k) = J4(i, k) + V(n2) * ...
                              (G(m,n2) * sin(del(m) - del(n2)) - ...
                               B(m,n2) * cos(del(m) - del(n2)));
                end
                J4(i, k) = J4(i, k) - V(m) * B(m, m);
            else
                J4(i, k) = V(m) * ...
                          (G(m,n) * sin(del(m) - del(n)) - ...
                           B(m,n) * cos(del(m) - del(n)));
            end
        end
    end
    
    % Assemble full Jacobian
    J = [J1 J2; J3 J4];
end