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
