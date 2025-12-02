% Experiment 6: Bad Data Detection and Elimination
clear
clc
close all


%% ini data
V = [1.06; 1.045; 1.01; 1.01772371530645; 1.01956929175548; 1.07; ...
     1.06158474200345; 1.09; 1.05599579866432; 1.05104189972113; ...
     1.05694910281384; 1.05520754822189; 1.05040466613156; 1.03557823098035];
delta = zeros(14,1);
B = [0.0528 0.0492 0.0438 0.0340 0.0346 0.0128 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
X = [0.05917; 0.22304; 0.19797; 0.17632; 0.17388; 0.17103; 0.04211; 0.1989; ...
     0.25581; 0.13027; 0.17615; 0.11001; 0.0845; 0.27038; 0.19207; 0.19988; ...
     0.34802; 0.20912; 0.55618; 0.20912];
R = [0.01938; 0.05403; 0.04699; 0.05811; 0.05695; 0.06701; 0.01335; 0.09498; ...
     0.12291; 0.06615; 0; 0; 0.03181; 0.12711; 0.08205; 0.22092; 0.17093; 0; 0; 0];
from_bus = [1; 1; 2; 2; 2; 3; 4; 4; 4; 5; 6; 6; 6; 7; 7; 9; 9; 10; 12; 13];
to_bus = [2; 5; 3; 4; 5; 4; 5; 7; 9; 6; 11; 12; 13; 8; 9; 10; 14; 11; 13; 14];

Z = R + 1j*X;
y = 1 ./ Z;
Nbus = length(V);
branches = length(from_bus);
[Ybuss] = Ybuss();
[z_true] = SE_Measurments();


rng(42); 
sigma = 0.02; 
noise = sigma * randn(length(z_true), 1);
z = z_true + noise;
z_original = z;

%% inject bad data
bad_measurement_index = 50; 
bad_data_magnitude = 0.15; 
z(bad_measurement_index) = z(bad_measurement_index) + bad_data_magnitude;

fprintf('Bad data  measurement index: %d\n', bad_measurement_index);
fprintf('Bad data : %.4f\n', bad_data_magnitude);
fprintf('Original value: %.6f, Corrupted value: %.6f\n\n', ...
    z_original(bad_measurement_index), z(bad_measurement_index));

W = eye(length(z)) / (sigma^2);

%% chi-squared 
alpha = 0.05;
m = length(z); 
n = 2*Nbus - 1; 
degrees_of_freedom = m - n;
chi2_threshold = chi2inv(1 - alpha, degrees_of_freedom);

fprintf('para:\n');
fprintf('  (alpha): %.2f\n', alpha);
fprintf('  m): %d\n', m);
fprintf('  n): %d\n', n);
fprintf('  Degrees of freedom: %d\n', degrees_of_freedom);
fprintf('   threshold: %.4f\n\n', chi2_threshold);

%% State Estimation with Bad Data Detection
iteration_outer = 0;
max_outer_iter = 10;
bad_data_detected = true;
bad_data_indices = [];
tau = 3.0; 
tau_adaptive = 0.15; 
active_measurements = true(length(z), 1); 

while bad_data_detected && iteration_outer < max_outer_iter
    iteration_outer = iteration_outer + 1;
    
    fprintf('----------------------------------------------------\n');
    fprintf('Outer Iteration %d: Running State Estimation\n', iteration_outer);
    

    x = [delta(2:end); V];
    z_active = z(active_measurements);
    W_active = W(active_measurements, active_measurements);
    % WLS 
    max_iter = 50;
    tolerance = 1e-4;
    iter = 0;
    converged = false;
    
    while iter < max_iter && ~converged
        iter = iter + 1;
        delta = [0; x(1:Nbus-1)];
        V = x(Nbus:end);
        [H_full] = State_Estimation_jacobian(V, delta, B, from_bus, to_bus, Ybuss);
        [f_full] = SEFunc(V, delta, y, from_bus, to_bus, Ybuss);
        H = H_full(active_measurements, :);
        f = f_full(active_measurements);
        e = z_active - f;
        G = (H' * W_active * H) \ (H' * W_active * e);%corr
        x = x + G;

        max_mismatch = max(abs(G));
        
        if max_mismatch < tolerance
            converged = true;
        end
    end
    
    fprintf('  State estimation converged in %d iterations\n', iter);
    
    %% Compute Final Residuals
    delta = [0; x(1:Nbus-1)];
    V = x(Nbus:end);
    [H_full] = State_Estimation_jacobian(V, delta, B, from_bus, to_bus, Ybuss);
    [f_full] = SEFunc(V, delta, y, from_bus, to_bus, Ybuss);
    
    % Sactive measurements
    H = H_full(active_measurements, :);
    f = f_full(active_measurements);
    r = z_active - f; % Residual vector
    
    %% Chi-Squared Test
    J_test = r' * W_active * r; 
    
    m_active = sum(active_measurements);
    degrees_of_freedom = m_active - n;
    chi2_threshold = chi2inv(1 - alpha, degrees_of_freedom);
    
    fprintf('  Chi-squared test statistic J(x): %.4f\n', J_test);
    fprintf('  Chi-squared threshold: %.4f\n', chi2_threshold);
    
    if J_test > chi2_threshold
        fprintf('  Result: BAD DATA SUSPECTED (J > threshold)\n\n');
        

        % = H(H'WH)^(-1)H'
        Omega = H * ((H' * W_active * H) \ H');
        
       % S = I - ΩW
        S = eye(m_active) - Omega * W_active;
        
        r_normalized = zeros(m_active, 1);
        for i = 1:m_active
            if S(i,i) > 1e-10 % Avoid division by zero
                r_normalized(i) = abs(r(i)) / sqrt(S(i,i));
            else
                r_normalized(i) = 0;
            end
        end
        
        [max_r_norm, bad_index_local] = max(r_normalized);
        
        % Map  
        active_indices = find(active_measurements);
        bad_index_original = active_indices(bad_index_local);
        
        %  top 5 
        [sorted_r_norm, sorted_idx] = sort(r_normalized, 'descend');
        fprintf('  Top 5 suspicious measurements:\n');
        for k = 1:min(5, length(sorted_r_norm))
            meas_idx = active_indices(sorted_idx(k));
            fprintf('    Measurement %3d: r_N = %.4f, residual = %.6f\n', ...
                meas_idx, sorted_r_norm(k), r(sorted_idx(k)));
        end
        fprintf('\n');
        
        fprintf('  Largest normalized residual: %.4f at measurement %d\n', max_r_norm, bad_index_original);
        
        % Adaptive stopping criterion based on Chi-squared test result
        J_ratio = J_test / chi2_threshold;
        
        fprintf('  Chi-squared ratio (J/threshold): %.2f\n', J_ratio);
        
        if max_r_norm > tau
            % Standard case: clear bad data with high normalized residual
            fprintf('  Action: REMOVING measurement %d (|r_N| = %.4f > %.1f)\n\n', bad_index_original, max_r_norm, tau);
            active_measurements(bad_index_original) = false;
            bad_data_indices = [bad_data_indices; bad_index_original];
            
        elseif J_ratio > 1.2 && max_r_norm > tau_adaptive
            fprintf('  Action: REMOVING measurement %d (J/threshold = %.2f, r_N = %.4f > %.2f)\n\n', ...
                bad_index_original, J_ratio, max_r_norm, tau_adaptive);
            active_measurements(bad_index_original) = false;
            bad_data_indices = [bad_data_indices; bad_index_original];
            
        else
            % No clear bad measurement
            if J_ratio > 1.1
                fprintf('  Warning: Chi-squared still elevated (%.2f × threshold) but no clear bad measurement\n', J_ratio);
                fprintf('  Action: Stopping to avoid removing good measurements\n\n');
            else
                fprintf('  Action: No bad measurement found. Stopping.\n\n');
            end
            bad_data_detected = false;
        end
    else
        fprintf('  Result: NO BAD DATA DETECTED (J <= threshold)\n\n');
        bad_data_detected = false;
    end
end

%% Final Results
fprintf('====================================================\n');
fprintf('   FINAL RESULTS (After Bad Data Removal)\n');
fprintf('====================================================\n\n');

fprintf('Bad measurements removed: %d\n', length(bad_data_indices));
if ~isempty(bad_data_indices)
    fprintf('Indices of removed measurements:\n');
    disp(bad_data_indices);
end

fprintf('\nFinal Estimated States:\n');
fprintf('Bus   Voltage (p.u.)   Angle (deg)\n');
fprintf('----------------------------------------\n');
delta_deg = rad2deg([0; x(1:Nbus-1)]);
for i = 1:Nbus
    fprintf('%2d    %8.6f       %8.4f\n', i, V(i), delta_deg(i));
end

%% Run comparison with clean data (no bad data injection)
fprintf('\n====================================================\n');
fprintf('   COMPARISON: Running with CLEAN data (no bad data)\n');
fprintf('====================================================\n\n');

% Reset measurements to clean version (with noise but no bad data)
z_clean = z_original;
active_clean = true(length(z_clean), 1);
x_clean = [zeros(Nbus-1,1); [1.06; 1.045; 1.01; 1.01772371530645; 1.01956929175548; 1.07; 1.06158474200345; 1.09; 1.05599579866432; 1.05104189972113; 1.05694910281384; 1.05520754822189; 1.05040466613156; 1.03557823098035]];

% Run one state estimation
max_iter = 50;
tolerance = 1e-4;
iter = 0;

while iter < max_iter
    iter = iter + 1;
    
    delta_clean = [0; x_clean(1:Nbus-1)];
    V_clean = x_clean(Nbus:end);
    
    [H_clean] = State_Estimation_jacobian(V_clean, delta_clean, B, from_bus, to_bus, Ybuss);
    [f_clean] = SEFunc(V_clean, delta_clean, y, from_bus, to_bus, Ybuss);
    
    e_clean = z_clean - f_clean;
    G_clean = (H_clean' * W * H_clean) \ (H_clean' * W * e_clean);
    x_clean = x_clean + G_clean;
    
    if max(abs(G_clean)) < tolerance
        break;
    end
end

% Compute Chi-squared for clean data
delta_clean = [0; x_clean(1:Nbus-1)];
V_clean = x_clean(Nbus:end);
[H_clean] = State_Estimation_jacobian(V_clean, delta_clean, B, from_bus, to_bus, Ybuss);
[f_clean] = SEFunc(V_clean, delta_clean, y, from_bus, to_bus, Ybuss);
r_clean = z_clean - f_clean;
J_clean = r_clean' * W * r_clean;

fprintf('Clean data Chi-squared statistic: %.4f (threshold: %.4f)\n', J_clean, chi2_threshold);
if J_clean <= chi2_threshold
    fprintf('Clean data result: PASS (no bad data detected)\n');
else
    fprintf('Clean data result: Still high (may need investigation)\n');
end
fprintf('Clean data converged in %d iterations\n\n', iter);

fprintf('Clean Data Estimated States:\n');
fprintf('Bus   Voltage (p.u.)   Angle (deg)\n');
fprintf('----------------------------------------\n');
delta_clean_deg = rad2deg([0; x_clean(1:Nbus-1)]);
for i = 1:Nbus
    fprintf('%2d    %8.6f       %8.4f\n', i, V_clean(i), delta_clean_deg(i));
end

fprintf('\n====================================================\n');
fprintf('   Algorithm completed successfully!\n');
fprintf('====================================================\n');