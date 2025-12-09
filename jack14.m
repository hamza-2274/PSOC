%% IEEE 14-Bus Power System State Estimation with Jacobian Calculation
clc; clear; close all;

%% ====================== SYSTEM DATA ======================
BMva = 100;
busdata = [
 1 1.060 0 232.4 -16.9 0.0 0.0 0 0 1;
 2 1.045 0 40.0 43.56 21.7 12.7 -40 50 2;
 3 1.010 0 0.0 0.0 94.2 19.0 -25 40 3;
 4 1.000 0 0.0 0.0 47.8 -3.9 0 0 3;
 5 1.000 0 0.0 0.0 7.6 1.6 0 0 3;
 6 1.070 0 0.0 0.0 11.2 7.5 -6 24 2;
 7 1.000 0 0.0 0.0 0.0 0.0 0 0 3;
 8 1.090 0 0.0 0.0 0.0 0.0 0 0 2;
 9 1.000 0 0.0 0.0 29.5 16.6 0 0 3;
 10 1.000 0 0.0 0.0 9.0 5.8 0 0 3;
 11 1.000 0 0.0 0.0 3.5 1.8 0 0 3;
 12 1.000 0 0.0 0.0 6.1 1.6 0 0 3;
 13 1.000 0 0.0 0.0 13.5 5.8 0 0 3;
 14 1.000 0 0.0 0.0 14.9 5.0 0 0 3;
];

linedata = [
 1 2 0.01938 0.05917 0.0264 1.000;
 1 5 0.05403 0.22304 0.0246 1.000;
 2 3 0.04699 0.19797 0.0219 1.000;
 2 4 0.05811 0.17632 0.0170 1.000;
 2 5 0.05695 0.17388 0.0173 1.000;
 3 4 0.06701 0.17103 0.0064 1.000;
 4 5 0.01335 0.04211 0.0000 1.000;
 4 7 0.00000 0.20912 0.0000 0.978;
 4 9 0.00000 0.55618 0.0000 0.969;
 5 6 0.00000 0.25202 0.0000 1.000;
 6 11 0.09498 0.19890 0.0000 1.000;
 6 12 0.12291 0.25581 0.0000 1.000;
 6 13 0.06615 0.13027 0.0000 1.000;
 7 8 0.00000 0.17615 0.0000 1.000;
 7 9 0.00000 0.11001 0.0000 1.000;
 9 10 0.03181 0.08450 0.0000 1.000;
 9 14 0.12711 0.27038 0.0000 1.000;
 10 11 0.08205 0.19207 0.0000 1.000;
 12 13 0.22092 0.19988 0.0000 1.000;
 13 14 0.17093 0.34802 0.0000 1.000;
];

nbus = size(busdata, 1);

%% ====================== Y-BUS FORMATION ======================
Ybuss = zeros(nbus, nbus);
Bsh = zeros(nbus, 1);  % Store shunt susceptances

for k = 1:size(linedata, 1)
    i = linedata(k, 1); 
    j = linedata(k, 2);
    R = linedata(k, 3); 
    X = linedata(k, 4); 
    Bc = linedata(k, 5); 
    tap = linedata(k, 6);
    
    if tap == 0, tap = 1; end
    
    Z = R + 1i * X; 
    if abs(Z) < 1e-10, Z = 1e-10; end
    
    y = 1 / Z; 
    ysh = 1i * Bc;
    
    Ybuss(i, j) = Ybuss(i, j) - y / tap; 
    Ybuss(j, i) = Ybuss(i, j);
    Ybuss(i, i) = Ybuss(i, i) + y / (tap^2) + ysh; 
    Ybuss(j, j) = Ybuss(j, j) + y + ysh;
    
    % Store total shunt susceptance for each line
    Bsh(k) = Bc;
end

%% ====================== LOAD OR INITIALIZE STATE ======================
if exist('V.txt', 'file') && exist('theta.txt', 'file')
    V = load('V.txt'); 
    delta = load('theta.txt');
else
    V = busdata(:, 2); 
    delta = busdata(:, 3) * pi / 180;
end

% Convert delta to radians if needed
if max(abs(delta)) > 2*pi
    delta = delta * pi / 180;
end

%% ====================== BRANCH DATA ======================
from_bus = linedata(:, 1);
to_bus = linedata(:, 2);
branches = length(from_bus);
%disp(to_bus);
% Create B vector (shunt susceptances)
B = Bsh;

%% ====================== ADMITTANCE PARAMETERS ======================
Nbus = length(V);
Ymag = abs(Ybuss);
Yang = angle(Ybuss);  % Keep in radians for calculations

%% ====================== JACOBIAN CALCULATIONS ======================

%% Jacobian J1
J1 = zeros(Nbus, Nbus-1);

%% Jacobian J2
J2 = eye(Nbus);

%% Jacobian J3: dQ/dδ
J3 = zeros(Nbus, Nbus-1);

for i = 1:Nbus
    for j = 2:Nbus
        if i ~= j  
            % Off-diagonal entries
            J3(i,j-1) = -V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
        else
            % Diagonal element
            sum_term = 0;
            for k = 1:Nbus
                if k ~= i
                    sum_term = sum_term + V(i)*V(k)*Ymag(i,k)*sin(Yang(i,k) + delta(k) - delta(i));
                end
            end
            J3(i,j-1) = sum_term;
        end
    end
end

%% Jacobian J4: dQ/dV
J4 = zeros(Nbus, Nbus);

for i = 1:Nbus
    for j = 1:Nbus
        if i ~= j  
            % Off-diagonal entries
            J4(i,j) = V(i)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i));
        else
            % Diagonal element
            diagonal_term = 2 * V(i) * V(i) * Ymag(i,i) * cos(Yang(i,i));
            sum_term = 0;
            for k = 1:Nbus
                if k ~= i
                    sum_term = sum_term + V(k)*Ymag(i,k)*cos(Yang(i,k) + delta(k) - delta(i));
                end
            end
            J4(i,j) = diagonal_term + sum_term; 
        end
    end
end

%% Jacobian J5: dP/dδ
J5 = zeros(Nbus, Nbus-1);

for i = 1:Nbus
    for j = 2:Nbus
        if i ~= j  
            % Off-diagonal entries
            J5(i,j-1) = -V(i)*V(j)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i));
        else
            % Diagonal element
            sum_term = 0;
            for k = 1:Nbus
                if k ~= i
                    sum_term = sum_term + V(i)*V(k)*Ymag(i,k)*cos(Yang(i,k) + delta(k) - delta(i));
                end
            end
            J5(i,j-1) = sum_term;
        end
    end
end

%% Jacobian J6: dP/dV
J6 = zeros(Nbus, Nbus);

for i = 1:Nbus
    for j = 1:Nbus
        if i ~= j  
            % Off-diagonal entries
            J6(i,j) = -V(i)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
        else
            % Diagonal element
            diagonal_term = -2 * V(i) * Ymag(i,i) * sin(Yang(i,i));
            sum_term = 0;
            for k = 1:Nbus
                if k ~= i
                    sum_term = sum_term - V(k)*Ymag(i,k)*sin(Yang(i,k) + delta(k) - delta(i));
                end
            end
            J6(i,j) = diagonal_term + sum_term; 
        end
    end
end

%% Jacobian J7: dPij/dδ (from end)
J7 = zeros(branches, Nbus-1);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    if i ~= 1
        J7(a,i-1) = V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
    end
    if j ~= 1
        J7(a,j-1) = -V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
    end
end

%% Jacobian J8: dPij/dV (from end)
J8 = zeros(branches, Nbus);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    J8(a,i) = (2*V(i)*real(Ybuss(i,j))) + (V(j)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i)));
    J8(a,j) = V(i)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i));
end

%% Jacobian J9: dPji/dδ (to end)
J9 = zeros(branches, Nbus-1);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    if i ~= 1
        J9(a,i-1) = -V(i)*V(j)*Ymag(j,i)*sin(Yang(i,j)+delta(i)-delta(j));
    end
    if j ~= 1
        J9(a,j-1) = V(i)*V(j)*Ymag(j,i)*sin(Yang(i,j)+delta(i)-delta(j));
    end
end

%% Jacobian J10: dPji/dV (to end)
J10 = zeros(branches, Nbus);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    J10(a,i) = V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(i)-delta(j));
    J10(a,j) = (2*V(j)*Ymag(j,i)*cos(Yang(j,i))) + (V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(i)-delta(j)));
end

%% Jacobian J11: dQij/dδ (from end)
J11 = zeros(branches, Nbus-1);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    if i ~= 1
        J11(a,i-1) = V(i)*V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(j)-delta(i));
    end
    if j ~= 1
        J11(a,j-1) = -V(i)*V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(j)-delta(i));
    end
end

%% Jacobian J12: dQij/dV (from end)
J12 = zeros(branches, Nbus);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    J12(a,j) = -V(i)*Ymag(j,i)*sin(Yang(j,i)+delta(j)-delta(i));
    BB = (B(a)/2) + imag(Ybuss(i,j));
    J12(a,i) = (-2*V(i)*BB) - (V(j)*Ymag(j,i)*sin(Yang(j,i)+delta(j)-delta(i)));
end

%% Jacobian J13: dQji/dδ (to end)
J13 = zeros(branches, Nbus-1);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    if i ~= 1
        J13(a,i-1) = -V(i)*V(j)*Ymag(j,i)*cos(Yang(j,i)+delta(i)-delta(j));
    end
    if j ~= 1
        J13(a,j-1) = V(i)*V(j)*Ymag(j,i)*cos(Yang(j,i)+delta(i)-delta(j));
    end
end

%% Jacobian J14: dQji/dV (to end)
J14 = zeros(branches, Nbus);

for a = 1:branches
    i = from_bus(a);
    j = to_bus(a);
    J14(a,i) = -V(j)*Ymag(j,i)*sin(Yang(j,i)+delta(i)-delta(j));
    BB = (B(a)/2) + imag(Ybuss(i,j));
    J14(a,j) = (-2*V(j)*BB) - (V(i)*Ymag(j,i)*sin(Yang(j,i)+delta(i)-delta(j)));
end

%% ====================== COMPLETE JACOBIAN MATRIX ======================
J = [J1  J2;
     J3  J4;
     J5  J6;
     J7  J8;
     J9  J10;
     J11 J12;
     J13 J14];

%% ====================== DISPLAY CALCULATED JACOBIANS ======================
fprintf('\n=============================================\n');
fprintf('CALCULATED JACOBIAN MATRICES\n');
fprintf('=============================================\n\n');

fprintf('J1 (Dimensions: %d x %d):\n', size(J1,1), size(J1,2));
disp(J1);

fprintf('\nJ2 (Dimensions: %d x %d):\n', size(J2,1), size(J2,2));
disp(J2);

fprintf('\nJ3 (Dimensions: %d x %d):\n', size(J3,1), size(J3,2));
disp(J3);

fprintf('\nJ4 (Dimensions: %d x %d):\n', size(J4,1), size(J4,2));
disp(J4);

fprintf('\nJ5 (Dimensions: %d x %d):\n', size(J5,1), size(J5,2));
disp(J5);

fprintf('\nJ6 (Dimensions: %d x %d):\n', size(J6,1), size(J6,2));
disp(J6);

fprintf('\nJ7 (Dimensions: %d x %d):\n', size(J7,1), size(J7,2));
disp(J7);

fprintf('\nJ8 (Dimensions: %d x %d):\n', size(J8,1), size(J8,2));
disp(J8);

fprintf('\nJ9 (Dimensions: %d x %d):\n', size(J9,1), size(J9,2));
disp(J9);

fprintf('\nJ10 (Dimensions: %d x %d):\n', size(J10,1), size(J10,2));
disp(J10);

fprintf('\nJ11 (Dimensions: %d x %d):\n', size(J11,1), size(J11,2));
disp(J11);

fprintf('\nJ12 (Dimensions: %d x %d):\n', size(J12,1), size(J12,2));
disp(J12);

fprintf('\nJ13 (Dimensions: %d x %d):\n', size(J13,1), size(J13,2));
disp(J13);

fprintf('\nJ14 (Dimensions: %d x %d):\n', size(J14,1), size(J14,2));
disp(J14);

%% ====================== ERROR STATISTICS ======================
fprintf('\n=============================================\n');
fprintf('ERROR STATISTICS\n');
fprintf('=============================================\n\n');

jacobian_names = {'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', ...
                  'J9', 'J10', 'J11', 'J12', 'J13', 'J14'};
jacobian_matrices = {J1, J2, J3, J4, J5, J6, J7, J8, J9, J10, J11, J12, J13, J14};

for k = 1:length(jacobian_names)
    jacobian_file = [jacobian_names{k} '.txt'];
    
    if exist(jacobian_file, 'file')
        J_ref = load(jacobian_file);
        J_calc = jacobian_matrices{k};
        
        if ~isequal(size(J_ref), size(J_calc))
            fprintf('%s: DIMENSION MISMATCH (Ref: %dx%d, Calc: %dx%d)\n\n', ...
                    jacobian_names{k}, size(J_ref,1), size(J_ref,2), ...
                    size(J_calc,1), size(J_calc,2));
            continue;
        end
        
        error_matrix = abs(J_ref - J_calc);
        
        fprintf('%s Error:\n', jacobian_names{k});
        fprintf('  Max Error:  %.6e\n', max(error_matrix(:)));
        fprintf('  Mean Error: %.6e\n', mean(error_matrix(:)));
        fprintf('  RMSE:       %.6e\n\n', sqrt(mean(error_matrix(:).^2)));
    else
        fprintf('%s: Reference file not found\n\n', jacobian_names{k});
    end
end

