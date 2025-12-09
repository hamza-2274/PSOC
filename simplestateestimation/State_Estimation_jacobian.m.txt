function [J] = State_Estimation_jacobian(V, delta,B,from_bus,to_bus,Ybuss)
[Ybuss]=Ybuss();

%V = [1.06; 1.045; 1.01; 1.01772371530645; 1.01956929175548; 1.07; 1.06158474200345; 1.09; 1.05599579866432; 1.05104189972113; 1.05694910281384; 1.05520754822189; 1.05040466613156; 1.03557823098035];
%delta = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
%B = [0.0528 0.0492 0.0438 0.0340 0.0346 0.0128 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
%from_bus = [1; 1; 2; 2; 2; 3; 4; 4; 4; 5; 6; 6; 6; 7; 7; 9; 9; 10; 12; 13];
%to_bus = [2; 5; 3; 4; 5; 4; 5; 7; 9; 6; 11; 12; 13; 8; 9; 10; 14; 11; 13; 14];





branches=length(from_bus); 
Nbus=length(V);
Ymag = abs(Ybuss);                 % magnitude
Yang = angle(Ybuss) ;      % angle in degrees
%% Jacobian no 1
J1=zeros(Nbus,Nbus-1);

%% Jacobian no 2
J2=eye(Nbus);


%% Jacobian no 3   \\done
J3 = zeros(Nbus,Nbus-1);

for i=1:Nbus
    for j=2:Nbus
            % off diagonal entries
        if i~=j  
           J3(i,j-1)= -V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
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

%% Jacobian no 4  \\done
J4 = zeros(Nbus,Nbus);

for i=1:Nbus
    for j=1:Nbus
            % off diagonal entries
        if i~=j  
           J4(i,j)= V(i)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i));
        else
             % Diagonal element
            diagonal_term = 2 * V(i)  * Ymag(i,i) * cos(Yang(i,i));  % NEW VARIABLE: diagonal_term
           
            sum_term = 0;
            for k = 1:Nbus
                if k ~= i
                   sum_term = sum_term + V(k)*Ymag(i,k)*cos(Yang(i,k) + delta(k)- delta(i));
                end
            end
            J4(i,j) = diagonal_term + sum_term; 
        end
    end
end

%% Jacobian no 5      \\done
J5 = zeros(Nbus,Nbus-1);

for i=1:Nbus
    for j=2:Nbus
            % off diagonal entries
        if i~=j  
           J5(i,j-1)= -V(i)*V(j)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i));
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

%% Jacobian no 6     \\ done
J6 = zeros(Nbus,Nbus);

for i=1:Nbus
    for j=1:Nbus
            % off diagonal entries
        if i~=j  
           J6(i,j)=- V(i)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
        else
             % Diagonal element
            diagonal_term = -2 * V(i)* Ymag(i,i) * sin(Yang(i,i));  % NEW VARIABLE: diagonal_term
           
            sum_term = 0;
            for k = 1:Nbus
                if k ~= i
                   sum_term = sum_term - V(k)*Ymag(i,k)*sin(Yang(i,k) + delta(k) - delta(i) );
                end
            end
            J6(i,j) = diagonal_term + sum_term; 
        end
    end
end

%% Jacobian no 7   
J7=zeros(branches,Nbus-1);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
      if i~=1
      J7(a,i-1)=V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
      end
      if j~=1
      J7(a,j-1)=-V(i)*V(j)*Ymag(i,j)*sin(Yang(i,j)+delta(j)-delta(i));
      end
end
%% Jacobian no 8
J8=zeros(branches,Nbus);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
    J8(a,i)=( -2*V(i)*real(Ybuss(i,j)) ) +  (V(j)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i)));
    J8(a,j)=V(i)*Ymag(i,j)*cos(Yang(i,j)+delta(j)-delta(i));
end

%% Jacobian no 9
J9=zeros(branches,Nbus-1);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
      if i~=1
      J9(a,i-1)=-V(i)*V(j)*Ymag(j,i)*sin(Yang(i,j)+delta(i)-delta(j));
      end
      if j~=1
      J9(a,j-1)=V(i)*V(j)*Ymag(j,i)*sin(Yang(i,j)+delta(i)-delta(j));
      end
end

%% Jacobian no 10
J10=zeros(branches,Nbus);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
    J10(a,i)=V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(i)-delta(j));
    J10(a,j)=(-2*V(j)*Ymag(j,i)*cos(Yang(j,i)) ) +(V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(i)-delta(j)));
end

%% Jacobian no 11
J11=zeros(branches,Nbus-1);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
      if i~=1
      J11(a,i-1)=V(i)*V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(j)-delta(i));
      end
      if j~=1
      J11(a,j-1)=-V(i)*V(j)*Ymag(j,i)*cos(Yang(i,j)+delta(j)-delta(i));
      end
end

%% Jacobian no 12    
J12=zeros(branches,Nbus);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
    J12(a,j)=-V(i)*Ymag(j,i)*sin(Yang(j,i)+delta(j)-delta(i));
    BB=(B(a)/2)- (imag( Ybuss (i,j) ) );
    J12(a,i)=(-2*V(i)*BB) - (V(j)*Ymag(j,i)*sin(Yang(j,i)+delta(j)-delta(i) ) );
    
end

%% Jacobian no 13
J13=zeros(branches,Nbus-1);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
      if i~=1
      J13(a,i-1)=-V(i)*V(j)*Ymag(j,i)*cos(Yang(j,i)+delta(i)-delta(j));
      end
      if j~=1
      J13(a,j-1)=V(i)*V(j)*Ymag(j,i)*cos(Yang(j,i)+delta(i)-delta(j));
      end
end

%% Jacobian no 14
J14=zeros(branches,Nbus);
for a=1:length(from_bus)
    i=from_bus(a);
    j=to_bus(a);
    J14(a,i)=-V(j)*Ymag(j,i)*sin(Yang(j,i)+delta(i)-delta(j));
    BB=(B(a)/2)- (imag( Ybuss (i,j) ) );
    J14(a,j)=(-2*V(j)*BB) - (V(i)*Ymag(j,i)*sin(Yang(j,i)+delta(i)-delta(j)));
    
end

J=[J1  J2;
   J3  J4;
   J5  J6;
   J7  J8;
   J9  J10;
   J11 J12;
   J13 J14];



end