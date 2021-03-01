function [x_der] = motion_eq(t,x,k,Vac,Vdc,w0,m,b,eps,A,g0)
    x_der(1,1) = x(2);
    x_der(2,1) = -k/m*x(1) - b/m*x_der(1,1) + ( eps*A*(Vac*cos(w0*t)+Vdc).^2 / (2*(g0-x(1))^2) )/m;
end