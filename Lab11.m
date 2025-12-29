clear all;
close all;
clc;

%% 1. System Parameters (Kundur Ex 12.7 / Manual)
H = 6.5;         
D = 20;          % Damping
Tdo_prime = 8.0; 
Tqo_prime = 0.4; 
xd = 1.8;       
xq = 1.7;       
xd_prime = 0.3;  
xq_prime = 0.55; 
Vt = 1.0;        % Terminal Voltage
freq = 60;  
w_syn = 2 * pi * freq; 

% Disturbance Parameters
t_disturb = 1.0; 
Pm0 = 0.8;       % Initial Mechanical Power
delta_Pm = 0.05 * Pm0; % 5% Step Increase

% Pack static parameters into structure
P = struct('H', H, 'D', D, 'Vt', Vt, ...
           'xd', xd, 'xq', xq, 'xd_prime', xd_prime, 'xq_prime', xq_prime, ...
           'Tdo_prime', Tdo_prime, 'Tqo_prime', Tqo_prime, ...
           'w_syn', w_syn, 'Pm0', Pm0, 't_disturb', t_disturb, 'delta_Pm', delta_Pm);

%% 2. Smart Initialization (Steady-State Calculation)
% We define the Operating Point first:
Pe0 = P.Pm0;      % Electrical Power must match Mechanical Power
Qe0 = 0.2;        % Assume a typical reactive power load (e.g. 0.2 pu)
Vt0 = P.Vt;       % Terminal Voltage

% Step A: Calculate Terminal Current Phasor (It)
% S = P + jQ = V * conj(I)  =>  I = conj(S/V)
S_complex = Pe0 + 1j*Qe0;
It_phasor = conj(S_complex / Vt0);

% Step B: Calculate Rotor Angle (delta)
% We use the quadrature axis reactance (xq) to find the internal angle
Eq_phasor = Vt0 + 1j * P.xq * It_phasor;
delta0 = angle(Eq_phasor); % This is the exact steady-state angle

% Step C: Park Transformation (d-q components) of Voltage & Current
% Rotate phasors by angle (delta - pi/2) to align with d-q frame
V_dq = Vt0 * exp(-1j * (delta0 - pi/2)); 
I_dq = It_phasor * exp(-1j * (delta0 - pi/2));

Vd0 = real(V_dq);
Vq0 = imag(V_dq);
Id0 = real(I_dq);
Iq0 = imag(I_dq);

% Step D: Calculate Steady-State Flux States (Eq', Ed')
Eq_prime0 = Vq0 + P.xd_prime * Id0; 
Ed_prime0 = Vd0 - P.xq_prime * Iq0; 

% Step E: Calculate Required Excitation (Ef)
Ef_required = Eq_prime0 + (P.xd - P.xd_prime) * Id0;

% UPDATE THE PARAMETER STRUCTURE
P.Ef = Ef_required; 
y0 = [delta0; 1.0; Eq_prime0; Ed_prime0]; 

fprintf('Initialization Results:\n');
fprintf('Required Excitation (Ef): %.4f pu\n', P.Ef);
fprintf('Initial Angle (delta): %.4f deg\n', delta0 * 180/pi);

%% 3. Simulation Setup
tspan = [0 10]; % Simulation time from 0 to 10 seconds
odefun = @(t, y) SMIB_f(t, y, P);

% Run Solver
[t, y] = ode45(odefun, tspan, y0);

%% 4. Plotting
figure('Color', 'w'); 

% Plot 1: Rotor Angle
subplot(3, 1, 1);
plot(t, y(:, 1) * 180 / pi, 'LineWidth', 1.5);
ylabel('Rotor Angle \delta (deg)');
title('SMIB Dynamic Response: 5% Step in P_m');
grid on;

% Plot 2: Speed Deviation
subplot(3, 1, 2);
plot(t, y(:, 2) - 1.0, 'r', 'LineWidth', 1.5); 
ylabel('Speed Dev \Delta\omega (pu)');
grid on;

% Plot 3: Internal Voltages
subplot(3, 1, 3);
plot(t, y(:, 3), 'b', 'LineWidth', 1.5, 'DisplayName', "E'_q");
hold on;
plot(t, y(:, 4), 'g--', 'LineWidth', 1.5, 'DisplayName', "E'_d");
ylabel('Internal Voltages (pu)');
xlabel('Time (s)');
legend('show', 'Location', 'best');
grid on;

%% 5. System Function
function dydt = SMIB_f(t, y, P)
    % Unpack State Variables
    delta = y(1);
    omega = y(2);
    Eq_prime = y(3);
    Ed_prime = y(4);
    
    % Apply Step Disturbance at t = 1.0s
    if t >= P.t_disturb
        Pm_current = P.Pm0 + P.delta_Pm;
    else
        Pm_current = P.Pm0;
    end
      
    % A. Algebraic Equations
    Vd = P.Vt * sin(delta);
    Vq = P.Vt * cos(delta);
    
    id = (Eq_prime - Vq) / P.xd_prime;
    iq = (Vd - Ed_prime) / P.xq_prime;
    
    % Electrical Power Output (Corrected for saliency)
    Pe = Vd * id + Vq * iq;
       
    % B. Differential Equations
    d_delta = P.w_syn * (omega - 1); 
    d_omega = (1 / (2 * P.H)) * (Pm_current - Pe - P.D * (omega - 1)); 
    d_Eq_prime = (1 / P.Tdo_prime) * (P.Ef - Eq_prime - (P.xd - P.xd_prime) * id); 
    d_Ed_prime = (1 / P.Tqo_prime) * (-Ed_prime + (P.xq - P.xq_prime) * iq); 
    
    dydt = [d_delta; d_omega; d_Eq_prime; d_Ed_prime];
end