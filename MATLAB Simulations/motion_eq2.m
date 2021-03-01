function [x_der] = motion_eq2(t,x,k,Vdc,m,b,eps,A,g0)
    x_der(1,1) = x(2);
    x_der(2,1) = -k/m*x(1) - b/m*x_der(1,1) + ( eps*A*Vdc^2 / (2*(g0-x(1))^2) )/m;
end