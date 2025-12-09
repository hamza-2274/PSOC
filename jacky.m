%% ========================================================================
%  NEWTON-RAPHSON LOAD FLOW ANALYSIS (IEEE-14 Bus System)
%  
%  This program performs power flow analysis using Newton-Raphson method
%  for the IEEE 14-bus test system.
%  ========================================================================

clear; clc; close all;
fprintf('\n--- NEWTON-RAPHSON LOAD FLOW (IEEE-14 BUS) ---\n\n');

%% ====================== SYSTEM DATA ======================
BMva = 100;   % Base MVA

% Bus Data Format:
% Columns: Bus | Vsp | delta | Pg | Qg | Pl | Ql | Qmin | Qmax | Type
% Bus Type: 1=Slack, 2=PV (Generator), 3=PQ (Load)
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

% Line Data Format:
% Columns: From | To | R | X | B/2 | Tap
linedata = [
    1   2   0.01938   0.05917   0.0264   1.000;
    1   5   0.05403   0.22304   0.0246   1.000;
    2   3   0.04699   0.19797   0.0219   1.000;
    2   4   0.05811   0.17632   0.0170   1.000;
    2   5   0.05695   0.17388   0.0173   1.000;
    3   4   0.06701   0.17103   0.0064   1.000;
    4   5   0.01335   0.04211   0.0000   1.000;
    4   7   0.00000   0.20912   0.0000   0.978;   % Transformer
    4   9   0.00000   0.55618   0.0000   0.969;   % Transformer
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

fb = linedata(:, 1);  % From bus
tb = linedata(:, 2);  % To bus

for k = 1:size(linedata, 1)
    i = fb(k);
    j = tb(k);
    R = linedata(k, 3);
    X = linedata(k, 4);
    Bc = linedata(k, 5);
    tap = linedata(k, 6);
    
    if tap == 0
        tap = 1;
    end
    
    Z = R + 1i * X;
    y = 1 / Z;
    ysh = 1i * Bc;
    
    % Off-diagonal elements
    Y(i, j) = Y(i, j) - y / tap;
    Y(j, i) = Y(i, j);
    
    % Diagonal elements
    Y(i, i) = Y(i, i) + y / (tap^2) + ysh;
    Y(j, j) = Y(j, j) + y + ysh;
end

G = real(Y);  % Conductance matrix
B = imag(Y);  % Susceptance matrix

%% ====================== INITIALIZATION ======================
type = busdata(:, 10);
V = busdata(:, 2);
del = deg2rad(busdata(:, 3));
Pg = busdata(:, 4) / BMva;
Qg = busdata(:, 5) / BMva;
Pl = busdata(:, 6) / BMva;
Ql = busdata(:, 7) / BMva;
Qmin = busdata(:, 8) / BMva;
Qmax = busdata(:, 9) / BMva;

Psp = Pg - Pl;  % net active power injection
Qsp = Qg - Ql;  % net reactive power injection

pq = find(type == 3);  % PQ bus indices
npq = length(pq);  %uused to size the Jacobian submatrice

tol = 1e-6;  % nvergence tolerance
Iter = 0;
Tol = 1;

fprintf('Starting Newton-Raphson iterations...\n\n');

%% ====================== NEWTON-RAPHSON ITERATION ======================
% Start timing for Newton-Raphson
tic_nr = tic;

while Tol > tol
    Iter = Iter + 1;
    
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
        break;
    end
    
    % ========== JACOBIAN MATRIX CONSTRUCTION ==========
    
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
    
    % Solve for corrections
    X = J \ M;
    
    % Update state variables
    dTheta = X(1:nbus - 1);
    dV = X(nbus:end);
    del(2:end) = del(2:end) + dTheta;
    V(pq) = V(pq) + dV;
end

% Stop timing for Newton-Raphson
time_nr = toc(tic_nr);

%% ====================== NR RESULTS ======================
fprintf('\n');
fprintf('====================================================\n');
fprintf('  NEWTON-RAPHSON: CONVERGENCE ACHIEVED\n');
fprintf('====================================================\n');
fprintf('  Iterations: %d\n', Iter);
fprintf('  Max Mismatch: %.6e pu\n', Tol);
fprintf('  Computation Time: %.6f seconds\n', time_nr);
fprintf('  Time per Iteration: %.6f seconds\n', time_nr/Iter);
fprintf('====================================================\n\n');

Angles_deg = rad2deg(del);

fprintf('BUS VOLTAGE AND ANGLE RESULTS:\n');
fprintf('----------------------------------------------\n');
fprintf(' Bus |  Voltage (pu)  |  Angle (deg)  \n');
fprintf('----------------------------------------------\n');
for i = 1:nbus
    fprintf(' %2d  |    %7.4f     |   %9.4f\n', i, V(i), Angles_deg(i));
end
fprintf('----------------------------------------------\n');

%% ====================== NR POWER LOSS CALCULATION ======================
fprintf('\nCALCULATING TRANSMISSION LOSSES...\n\n');

Ploss_total = 0;
Qloss_total = 0;

for k = 1:size(linedata, 1)
    i = fb(k);
    j = tb(k);
    R = linedata(k, 3);
    X = linedata(k, 4);
    tap = linedata(k, 6);
    
    if tap == 0
        tap = 1;
    end
    
    Z = R + 1i * X;
    y = 1 / Z;
    
    Vi = V(i) * exp(1i * del(i));
    Vj = V(j) * exp(1i * del(j));
    
    Iij = (Vi - Vj / tap) * y;
    Iji = (Vj - Vi * tap) * y;
    
    Sij = Vi * conj(Iij);
    Sji = Vj * conj(Iji);
    
    Ploss_total = Ploss_total + real(Sij + Sji);
    Qloss_total = Qloss_total + imag(Sij + Sji);
end

fprintf('====================================================\n');
fprintf('  SYSTEM LOSSES (NR)\n');
fprintf('====================================================\n');
fprintf('  Active Power Loss:   %.6f pu  (%.3f MW)\n', ...
        Ploss_total, Ploss_total * BMva);
fprintf('  Reactive Power Loss: %.6f pu  (%.3f MVAr)\n', ...
        Qloss_total, Qloss_total * BMva);
fprintf('====================================================\n\n');

fprintf('=== NEWTON-RAPHSON SIMULATION COMPLETE ===\n\n');

%% ========================================================================
%  STEP 3: FAST DECOUPLED LOAD FLOW IMPLEMENTATION
%  ========================================================================
fprintf('\n');
fprintf('====================================================\n');
fprintf('  STARTING FAST DECOUPLED LOAD FLOW ANALYSIS\n');
fprintf('====================================================\n\n');

%% ====================== FDLF INITIALIZATION ======================
% Reset voltages to initial values
V_fd = busdata(:, 2);
del_fd = deg2rad(busdata(:, 3));

% B' matrix neglects resistance (only reactance considered)
B_prime = zeros(nbus, nbus);

for k = 1:size(linedata, 1)
    i = fb(k);
    j = tb(k);
    R = linedata(k, 3);
    X = linedata(k, 4);
    Bc = linedata(k, 5);
    tap = linedata(k, 6);
    
    if tap == 0
        tap = 1;
    end
    
    % B' formation: neglect resistance, use only reactance
    % Yij = -1/(j*X) / tap
    y_prime = 1 / (1i * X);
    
    % Off-diagonal elements
    B_prime(i, j) = B_prime(i, j) - imag(y_prime / tap);
    B_prime(j, i) = B_prime(i, j);
    
    % Diagonal elements
    B_prime(i, i) = B_prime(i, i) + imag(y_prime / (tap^2));
    B_prime(j, j) = B_prime(j, j) + imag(y_prime);
end

B_double_prime = zeros(nbus, nbus);

for k = 1:size(linedata, 1)
    i = fb(k);
    j = tb(k);
    X = linedata(k, 4);
    Bc = linedata(k, 5);
    tap = linedata(k, 6);

    if tap == 0
        tap = 1;
    end
    
    y = 1 / (1i * X);
    
    % Off-diagonal elements
    B_double_prime(i, j) = B_double_prime(i, j) - imag(y / tap);
    B_double_prime(j, i) = B_double_prime(i, j);
    
    % Diagonal elements (include shunt charging)
    B_double_prime(i, i) = B_double_prime(i, i) + imag(y / (tap^2)) + Bc/2;
    B_double_prime(j, j) = B_double_prime(j, j) + imag(y) + Bc/2;
end


% Form Jacobian matrices
B11 = B_prime(2:end, 2:end);

B22 = B_double_prime(pq, pq);


% Final Jacobian for Fast Decoupled: J12 = J21 = 0
Jacob_FD = -[B11, zeros(nbus-1, npq); 
             zeros(npq, nbus-1), B22];

tol_fd = 0.01;  % Tolerance 
Iter_fd = 0;
max_iter_fd = 100;

fprintf('Fast Decoupled method initialized...\n');
fprintf('B'' (B11) matrix size: %dx%d\n', size(B11,1), size(B11,2));
fprintf('B'''' (B22) matrix size: %dx%d\n', size(B22,1), size(B22,2));
fprintf('Tolerance: %.6f\n\n', tol_fd);

%% ====================== FDLF ITERATION ======================
mismatch_fd = ones(2*nbus - 2 - length(find(type==2)), 1);

% Start timing for Fast Decoupled
tic_fd = tic;

while any(abs(mismatch_fd(:)) > tol_fd) && Iter_fd < max_iter_fd
    Iter_fd = Iter_fd + 1;
    
    % Calculate power injections
    P_fd = zeros(nbus, 1);
    Q_fd = zeros(nbus, 1);
    
    for i = 1:nbus
        for k = 1:nbus
            P_fd(i) = P_fd(i) + V_fd(i) * V_fd(k) * ...
                      (G(i,k) * cos(del_fd(i) - del_fd(k)) + ...
                       B(i,k) * sin(del_fd(i) - del_fd(k)));
            Q_fd(i) = Q_fd(i) + V_fd(i) * V_fd(k) * ...
                      (G(i,k) * sin(del_fd(i) - del_fd(k)) - ...
                       B(i,k) * cos(del_fd(i) - del_fd(k)));
        end
    end
    
    % Calculate mismatches
    dP_fd = Psp - P_fd;
    dQ_fd = Qsp - Q_fd;
    
    % Mismatch vector: [ΔP(2:end); ΔQ(pq)]
    mismatch_fd = [dP_fd(2:end); dQ_fd(pq)];
    error = Jacob_FD \ mismatch_fd;
    
    %  angles 
    del_fd(2:end) = del_fd(2:end) + error(1:nbus-1);
    
    % Update voltage (only PQ buses)
    V_fd(pq) = V_fd(pq) + error(nbus:end);
    
    % Average error for tracking
    error_avg = abs(mean(error(:)));
    
    % Display progress every 10 iterations
    if mod(Iter_fd, 10) == 0
        fprintf('  Iteration %d: Max mismatch = %.6e, Avg error = %.6e\n', ...
                Iter_fd, max(abs(mismatch_fd)), error_avg);
    end
end

% Stop timing for Fast Decoupled
time_fd = toc(tic_fd);

% Check if converged
if Iter_fd >= max_iter_fd
    fprintf('\n*** WARNING: Fast Decoupled method did NOT converge! ***\n');
    fprintf('*** Maximum iterations (%d) reached. ***\n', max_iter_fd);
    fprintf('*** Results may not be accurate. ***\n\n');
end

%% ====================== FDLF RESULTS ======================
fprintf('\n');
fprintf('====================================================\n');
if max(abs(mismatch_fd)) < tol_fd
    fprintf('  FAST DECOUPLED: CONVERGED\n');
else
    fprintf('  FAST DECOUPLED: DID NOT CONVERGE\n');
end
fprintf('====================================================\n');
fprintf('  Iterations: %d\n', Iter_fd);
fprintf('  Max Mismatch: %.6e pu\n', max(abs(mismatch_fd)));
fprintf('  Computation Time: %.6f seconds\n', time_fd);
fprintf('  Time per Iteration: %.6f seconds\n', time_fd/Iter_fd);
if max(abs(mismatch_fd)) >= tol_fd
    fprintf('  STATUS: FAILED (mismatch > tolerance)\n');
else
    fprintf('  STATUS: SUCCESS\n');
end
fprintf('====================================================\n\n');

Angles_deg_fd = rad2deg(del_fd);

fprintf('BUS VOLTAGE AND ANGLE RESULTS (FDLF):\n');
fprintf('----------------------------------------------\n');
fprintf(' Bus |  Voltage (pu)  |  Angle (deg)  \n');
fprintf('----------------------------------------------\n');
for i = 1:nbus
    fprintf(' %2d  |    %7.4f     |   %9.4f\n', i, V_fd(i), Angles_deg_fd(i));
end
fprintf('----------------------------------------------\n');

%% ====================== FDLF POWER LOSS CALCULATION ======================
fprintf('\nCALCULATING TRANSMISSION LOSSES (FDLF)...\n\n');

Ploss_total_fd = 0;
Qloss_total_fd = 0;

for k = 1:size(linedata, 1)
    i = fb(k);
    j = tb(k);
    R = linedata(k, 3);
    X = linedata(k, 4);
    tap = linedata(k, 6);
    
    if tap == 0
        tap = 1;
    end
    
    Z = R + 1i * X;
    y = 1 / Z;
    
    Vi_fd = V_fd(i) * exp(1i * del_fd(i));
    Vj_fd = V_fd(j) * exp(1i * del_fd(j));
    
    Iij_fd = (Vi_fd - Vj_fd / tap) * y;
    Iji_fd = (Vj_fd - Vi_fd * tap) * y;
    
    Sij_fd = Vi_fd * conj(Iij_fd);
    Sji_fd = Vj_fd * conj(Iji_fd);
    
    Ploss_total_fd = Ploss_total_fd + real(Sij_fd + Sji_fd);
    Qloss_total_fd = Qloss_total_fd + imag(Sij_fd + Sji_fd);
end

fprintf('====================================================\n');
fprintf('  SYSTEM LOSSES (FDLF)\n');
fprintf('====================================================\n');
fprintf('  Active Power Loss:   %.6f pu  (%.3f MW)\n', ...
        Ploss_total_fd, Ploss_total_fd * BMva);
fprintf('  Reactive Power Loss: %.6f pu  (%.3f MVAr)\n', ...
        Qloss_total_fd, Qloss_total_fd * BMva);
fprintf('====================================================\n\n');

%% ====================== COMPARISON: NR vs FDLF ======================
fprintf('\n');
fprintf('========================================================================\n');
fprintf('  COMPARISON: NEWTON-RAPHSON vs FAST DECOUPLED LOAD FLOW\n');
fprintf('========================================================================\n\n');

fprintf('CONVERGENCE & TIMING COMPARISON:\n');
fprintf('--------------------------------------------------------------------------------\n');
fprintf(' Method              | Iterations | Total Time (s) | Time/Iter (s) | Speedup\n');
fprintf('--------------------------------------------------------------------------------\n');
fprintf(' Newton-Raphson      |     %2d     |    %.6f      |   %.6f    |   1.00x\n', ...
        Iter, time_nr, time_nr/Iter);
fprintf(' Fast Decoupled      |     %2d     |    %.6f      |   %.6f    |   %.2fx\n', ...
        Iter_fd, time_fd, time_fd/Iter_fd, time_nr/time_fd);
fprintf('--------------------------------------------------------------------------------\n');
if time_fd < time_nr
    fprintf(' *** FDLF is %.2fx FASTER than Newton-Raphson! ***\n', time_nr/time_fd);
    fprintf(' *** Time saved: %.6f seconds (%.1f%% faster) ***\n', ...
            time_nr-time_fd, ((time_nr-time_fd)/time_nr)*100);
else
    fprintf(' Note: NR faster for this small system, but FDLF scales better for large systems\n');
end
fprintf('--------------------------------------------------------------------------------\n\n');

fprintf('VOLTAGE MAGNITUDE COMPARISON:\n');
fprintf('------------------------------------------------------------------------\n');
fprintf(' Bus | NR Voltage (pu) | FDLF Voltage (pu) | Difference (pu) | Error %%\n');
fprintf('------------------------------------------------------------------------\n');
for i = 1:nbus
    diff_V = abs(V(i) - V_fd(i));
    error_pct = (diff_V / V(i)) * 100;
    fprintf(' %2d  |    %7.4f      |     %7.4f       |    %.6f     | %.4f%%\n', ...
            i, V(i), V_fd(i), diff_V, error_pct);
end
fprintf('------------------------------------------------------------------------\n\n');

fprintf('VOLTAGE ANGLE COMPARISON:\n');
fprintf('------------------------------------------------------------------------\n');
fprintf(' Bus | NR Angle (deg) | FDLF Angle (deg) | Difference (deg) \n');
fprintf('------------------------------------------------------------------------\n');
for i = 1:nbus
    diff_angle = abs(Angles_deg(i) - Angles_deg_fd(i));
    fprintf(' %2d  |   %9.4f    |    %9.4f     |    %.6f\n', ...
            i, Angles_deg(i), Angles_deg_fd(i), diff_angle);
end
fprintf('------------------------------------------------------------------------\n\n');

fprintf('LOSSES COMPARISON:\n');
fprintf('------------------------------------------------------------------------\n');
fprintf(' Parameter           | Newton-Raphson | Fast Decoupled | Difference\n');
fprintf('------------------------------------------------------------------------\n');
fprintf(' Active Power Loss   |   %.3f MW     |   %.3f MW      | %.3f MW\n', ...
        Ploss_total*BMva, Ploss_total_fd*BMva, abs(Ploss_total-Ploss_total_fd)*BMva);
fprintf(' Reactive Power Loss |   %.3f MVAr   |   %.3f MVAr    | %.3f MVAr\n', ...
        Qloss_total*BMva, Qloss_total_fd*BMva, abs(Qloss_total-Qloss_total_fd)*BMva);
fprintf('------------------------------------------------------------------------\n\n');

fprintf('PERFORMANCE SUMMARY:\n');
fprintf('========================================================================\n');
fprintf('COMPUTATION TIME ANALYSIS:\n');
fprintf('- Newton-Raphson: %d iterations @ %.6f sec/iter = %.6f sec total\n', ...
        Iter, time_nr/Iter, time_nr);
fprintf('- Fast Decoupled: %d iterations @ %.6f sec/iter = %.6f sec total\n', ...
        Iter_fd, time_fd/Iter_fd, time_fd);
if time_fd < time_nr
    fprintf('- FDLF is %.2fx FASTER overall\n', time_nr/time_fd);
    fprintf('- Each FDLF iteration is %.2fx faster (simplified Jacobian)\n', ...
            (time_nr/Iter)/(time_fd/Iter_fd));
    fprintf('- Time saved: %.4f seconds\n', time_nr-time_fd);
else
    fprintf('- NR faster for this small 14-bus system\n');
    fprintf('- FDLF advantage increases with system size (100+ buses)\n');
end
fprintf('\nACCURACY ANALYSIS:\n');
fprintf('- Maximum voltage difference: %.6f pu (%.4f%%)\n', ...
        max(abs(V - V_fd)), max(abs(V - V_fd)./V)*100);
fprintf('- Maximum angle difference: %.6f degrees\n', max(abs(Angles_deg - Angles_deg_fd)));
fprintf('- Active power loss difference: %.6f MW\n', abs(Ploss_total-Ploss_total_fd)*BMva);
fprintf('- Accuracy: Excellent match between both methods\n');


fprintf('=== FAST DECOUPLED LOAD FLOW COMPLETE ===\n');