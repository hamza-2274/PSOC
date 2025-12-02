function [b] = SE_Measurments()
V = [1.06; 1.045; 1.01; 1.01772371530645; 1.01956929175548; 1.07; 1.06158474200345; 1.09; 1.05599579866432; 1.05104189972113; 1.05694910281384; 1.05520754822189; 1.05040466613156; 1.03557823098035];
delta = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
X = [0.05917; 0.22304; 0.19797; 0.17632; 0.17388; 0.17103; 0.04211; 0.1989; 0.25581; 0.13027; 0.17615; 0.11001; 0.0845; 0.27038; 0.19207; 0.19988; 0.34802; 0.20912; 0.55618; 0.20912];
R = [0.01938; 0.05403; 0.04699; 0.05811; 0.05695; 0.06701; 0.01335; 0.09498; 0.12291; 0.06615; 0; 0; 0.03181; 0.12711; 0.08205; 0.22092; 0.17093; 0; 0; 0];
B = [0.0528 0.0492 0.0438 0.0340 0.0346 0.0128 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
from_bus = [1; 1; 2; 2; 2; 3; 4; 4; 4; 5; 6; 6; 6; 7; 7; 9; 9; 10; 12; 13];
to_bus = [2; 5; 3; 4; 5; 4; 5; 7; 9; 6; 11; 12; 13; 8; 9; 10; 14; 11; 13; 14];
bus_type = [3; 2; 2; 1; 1; 2; 1; 2; 1; 1; 1; 1; 1; 1];
PGen = [0; 0.4; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
QGen = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
QLoad = [0; 0.127; 0.19; -0.039; 0.016; 0.075; 0; 0; 0.166; 0.058; 0.018; 0.016; 0.058; 0.05];
PLoad = [0; 0.217; 0.942; 0.478; 0.076; 0.112; 0; 0; 0.295; 0.09; 0.035; 0.061; 0.135; 0.149];
Z = R + 1j*X;     % complex impedance vector
y = 1 ./ Z;       % element-wise reciprocal to get admittance

[Ybus]=Ybuss();
Ymag = abs(Ybus);                 % magnitude
Yang = angle(Ybus);               % angle in radians
branches=length(from_bus);
Nbus = length(V);               % total no of buses
nPV = sum(bus_type == 2);       % PV buses
nPQ = sum(bus_type == 1);       % PQ buses
% Define PV and PQ buses
pv = find(bus_type == 2);
pq = find(bus_type == 1);
% Create new variable: list of all non-slack buses
non_slack_buses = nPV+nPQ;

PN=PGen-PLoad;
QN=QGen-QLoad;

% Iteration parameters
max_iter = 50;
tolerance = 1e-4;
iter = 0;
converged = false;
while iter < max_iter && ~converged
    iter = iter + 1;

    % Compute P injections
    P = zeros(Nbus,1);
    for i = 1:Nbus
        for j = 1:Nbus
            P(i) = P(i) + V(i)*V(j)*Ymag(i,j)*cos(Yang(i,j) - delta(i) + delta(j));
        end
    end
    
    % Compute Q injections
    Q = zeros(Nbus,1);
    for i = 1:Nbus
        for j = 1:Nbus
            Q(i) = Q(i) - V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j) - delta(i) + delta(j));
        end
    end

    M1=PN-P;
    M2=QN-Q;
    MM=[M1(2:end);M2(pq)];
    
    % Check convergence
    max_mismatch = max(abs(MM));
    fprintf('%4d      %12.6f\n', iter, max_mismatch);  
   
    if max_mismatch < tolerance
        converged = true;
        break;
    end

    J1 = zeros((nPQ+nPV), (nPQ+nPV));
    for i=2:Nbus
        for j=2:Nbus
            if i~=j  
               J1(i-1,j-1)= -V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)-delta(i)+delta(j));
            else
                sum_term = 0;
                for k = 1:Nbus
                    if k ~= i
                        sum_term = sum_term + V(i)*V(k)*Ymag(i,k)*sin(Yang(i,k) - delta(i) + delta(k));
                    end
                end
                J1(i-1,i-1) = sum_term;
            end
        end
    end

    J2 = zeros((nPQ+nPV), (nPQ));
    for a = 1:length(non_slack_buses)
        i = non_slack_buses(a);
        for b = 1:length(pq)
            j = pq(b);
            if i ~= j
                J2(a,b) = V(i)*Ymag(i,j)*cos(Yang(i,j) - delta(i) + delta(j));
            else
                diagonal_term = 2 * V(i) * Ymag(i,i) * cos(Yang(i,i));
                sum_term = 0;
                for k = 1:Nbus
                    if k ~= i  
                        sum_term = sum_term + V(k)*Ymag(i,k)*cos(Yang(i,k) - delta(i) + delta(k));
                    end
                end
                J2(a,b) = diagonal_term + sum_term; 
            end
        end
    end

    J3 = zeros(nPQ, (nPQ + nPV));
    for a = 1:length(pq)
        i = pq(a);
        for b = 1:length(non_slack_buses)
            j = non_slack_buses(b);
            if i ~= j
                J3(a,b) = - V(i) * V(j) * Ymag(i,j) *cos(Yang(i,j) - delta(i) + delta(j));
            else
                sum_term = 0;
                for k = 1:Nbus
                    if k ~= i
                        sum_term = sum_term + V(i)*V(k)*Ymag(i,k) * cos(Yang(i,k) - delta(i) + delta(k));
                    end
                end
                J3(a,b) = sum_term;
            end
        end
    end

    J4 = zeros(nPQ, nPQ);
    for a = 1:length(pq)
        i = pq(a);
        for b = 1:length(pq)
            j = pq(b);
            if i ~= j
                J4(a,b) = -V(i) * Ymag(i,j) * sin(Yang(i,j) - delta(i) + delta(j));
            else
                sum_term = 0;
                for k = 1:Nbus
                    if k ~= i
                        sum_term = sum_term + V(k)*Ymag(i,k)*sin(Yang(i,k) - delta(i) + delta(k));
                    end
                end
                J4(a,b) = -2*V(i)*Ymag(i,i)*sin(Yang(i,i)) - sum_term;
            end
        end
    end

    J = [J1  J2;
         J3  J4];

    Q_update = J\MM;
    delta(2:end) = delta(2:end) + Q_update(1:(Nbus-1));
    V(pq) = V(pq) + Q_update(Nbus:end);
end

% Final calculation of P and Q (CORRECTED - no adding to existing values)
P = zeros(Nbus,1);
Q = zeros(Nbus,1);
for i = 1:Nbus
    for j = 1:Nbus
        P(i) = P(i) + V(i)*V(j)*Ymag(i,j)*cos(Yang(i,j) - delta(i) + delta(j));
        Q(i) = Q(i) - V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j) - delta(i) + delta(j));
    end
end

% Branch Power Flow
V_complex = V .* exp(1j * delta);
S=zeros(branches,1);
S1=zeros(branches,1);
for a=1:branches
    i=from_bus(a);
    j=to_bus(a);
    S(a)=V_complex(i)*conj( y(a)* ( V_complex(i)-V_complex(j) ) );
    S1(a)=V_complex(j)*conj( y(a)* ( V_complex(j)-V_complex(i) ) );
end
P1=real(S);
Q1=imag(S);
P2=real(S1);
Q2=imag(S1);

% Assemble measurement vector
b=zeros(3*Nbus+4*branches,1);
b(1:Nbus)=V;                                    % Voltage magnitudes
b(Nbus+1:2*Nbus)=P;                             % Power injections
b(2*Nbus+1:3*Nbus)=Q;                           % Reactive injections
b(3*Nbus+1:3*Nbus+branches)=P1;                 % P from-bus
b(3*Nbus+branches+1:3*Nbus+2*branches)=P2;      % P to-bus
b(3*Nbus+2*branches+1:3*Nbus+3*branches)=Q1;    % Q from-bus
b(3*Nbus+3*branches+1:3*Nbus+4*branches)=Q2;    % Q to-bus

end