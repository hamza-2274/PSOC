function [Ybuss] = Ybuss()
    % Y-BUS FORMATION FUNCTION
    % Returns the admittance matrix and its real and imaginary components
    
    % Base MVA
    BMva = 100;
    
    % Line Data Format:
    % Columns: from_bus | to_bus | R | X | B/2 | tap
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
    
    % Determine number of buses
    nbus = max(max(linedata(:, 1:2)));
    
    % Initialize Y-bus matrix
    Ybuss = zeros(nbus, nbus);
    
    % Extract line data
    fb = linedata(:, 1);  % from bus
    tb = linedata(:, 2);  % to bus
    
    % Build Y-bus
    for k = 1:size(linedata, 1)
        i = fb(k);
        j = tb(k);
        R = linedata(k, 3);
        X = linedata(k, 4);
        Bc = linedata(k, 5);
        tap = linedata(k, 6);
        
        % Handle transformer tap (if tap = 0, set to 1)
        if tap == 0
            tap = 1;
        end
        
        % Calculate series impedance and admittance
        Z = R + 1i * X;
        y = 1 / Z;
        
        % Shunt admittance (line charging)
        ysh = 1i * Bc;
        
        % Off-diagonal elements (mutual admittance)
        Ybuss(i, j) = Ybuss(i, j) - y / tap;
        Ybuss(j, i) = Ybuss(i, j);
        
        % Diagonal elements (self admittance)
        Ybuss(i, i) = Ybuss(i, i) + y / (tap^2) + ysh;
        Ybuss(j, j) = Ybuss(j, j) + y + ysh;
    end
    
    % Extract real and imaginary parts
    G = real(Ybuss);
    B = imag(Ybuss);
end