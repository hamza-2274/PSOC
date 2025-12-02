clear
clc
tic
V = [1.06; 1.045; 1.01; 1.01772371530645; 1.01956929175548; 1.07; 1.06158474200345; 1.09; 1.05599579866432; 1.05104189972113; 1.05694910281384; 1.05520754822189; 1.05040466613156; 1.03557823098035];
delta = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
B = [0.0528 0.0492 0.0438 0.0340 0.0346 0.0128 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
X = [0.05917; 0.22304; 0.19797; 0.17632; 0.17388; 0.17103; 0.04211; 0.1989; 0.25581; 0.13027; 0.17615; 0.11001; 0.0845; 0.27038; 0.19207; 0.19988; 0.34802; 0.20912; 0.55618; 0.20912];
R = [0.01938; 0.05403; 0.04699; 0.05811; 0.05695; 0.06701; 0.01335; 0.09498; 0.12291; 0.06615; 0; 0; 0.03181; 0.12711; 0.08205; 0.22092; 0.17093; 0; 0; 0];
from_bus = [1; 1; 2; 2; 2; 3; 4; 4; 4; 5; 6; 6; 6; 7; 7; 9; 9; 10; 12; 13];
to_bus = [2; 5; 3; 4; 5; 4; 5; 7; 9; 6; 11; 12; 13; 8; 9; 10; 14; 11; 13; 14];
Z = R + 1j*X;     % complex impedance vector
y = 1 ./ Z;       % element-wise reciprocal to get admittance

Nbus=length(V);

x=[delta(2:end);V];         % state vector
[Ybuss]=Ybuss();
[b1] = SE_Measurments();     % measurment vector
%disp(H);
%disp(f); 
r = 0 + 0.02.*randn(122,1);
b=b1+r;
W=eye(122)/ (0.02^2);

% Iteration parameters
max_iter = 70;
tolerance = 1e-4;
iter = 0;
converged = false;
while iter < max_iter && ~converged
    iter = iter + 1;
    %ite7
[H] =  State_Estimation_jacobian(V, delta,B,from_bus,to_bus,Ybuss);
[f] =  SEFunc(V, delta,y,from_bus,to_bus,Ybuss);

e=b-f;

 G = (H' * W * H) \ (H' * W * e);
 x = x + G;

V=[x(Nbus:end)];
delta=[0;x(1:Nbus-1)];
% Check convergence
    max_mismatch = max(abs(G));
   
   fprintf('%4d      %12.6f\n', iter, max_mismatch);  
   
    if max_mismatch < tolerance
        converged = true;
        break      
    end
end
delta1=rad2deg(x(1:Nbus-1))

toc




