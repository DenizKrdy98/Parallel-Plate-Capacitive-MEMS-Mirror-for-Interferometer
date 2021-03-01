% *Program Initialization*
clear all;
close all;
clc

%% Input Variables
g0 = 300e-6;            % initial gap
gmin = 60e-6;           % minimum gap
W = 50e-3;              % width of plate
L = 50e-3;              % length of plate
T = 0.1e-6;             % thickness of plate
f0 = 100;               % resonant frequency (Hz)
Q = 0.4;                  % Q factor
d_plate = 2330;         % density of plate
eps = 8.854e-12;

% Mirror Parameters
d_mir = 2700;           % density of mirror
Wmir = 1e-3;            % width of mirror
Lmir = 1e-3;            % length of mirror
Tmir = 10e-6;          % thickness of mirror
m_mir = Wmir*Lmir*Tmir*d_mir;   % mass of mirror


%% Dependent Variables
w0 = 2*pi*f0;                     % rad/sec
a = w0/(2*Q);                     % damping coef
A = W*L;                          % Area of plate
m = d_plate*(A*T) + m_mir;        % mass of plate + mirror
k = m * w0^2;                     % spring coef
b = 2*m*a;                        % damper coef

%% ODE related Values
t_f = 400e-3;
x_i = 0;
xmax = 30e-6; 
tol = 1e-24;
tol_option = odeset('AbsTol',tol);

%% ODE Solution
Vdc1 = 7;
Vdc2 = 9;
Vdc3 = 9.3;
Vpi = 9.34;
Vdc4 = 9.7;

[time1,xt] = ode45(@(t,x) motion_eq2(t,x,k,Vdc1,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
x1 = xt(:,1);
[time2,xt] = ode45(@(t,x) motion_eq2(t,x,k,Vdc2,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
x2 = xt(:,1);
[time3,xt] = ode45(@(t,x) motion_eq2(t,x,k,Vdc3,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
x3 = xt(:,1);
[timepi,xt] = ode45(@(t,x) motion_eq2(t,x,k,Vpi,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
xpi = xt(:,1);
[time4,xt] = ode45(@(t,x) motion_eq2(t,x,k,Vdc4,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
x4 = xt(:,1);

for i = 1:length(x3)
    if(x3(i)>g0-gmin)
        x3(i)=g0-gmin;
    end
end
for i = 1:length(xpi)
    if(xpi(i)>g0-gmin)
        xpi(i)=g0-gmin;
    end
end
for i = 1:length(x4)
    if(x4(i)>g0-gmin)
        x4(i)=g0-gmin;
    end
end

figure
p1 = plot(time1,x1*10^6,'LineWidth',1);
hold on 
p2 = plot(time2,x2*10^6,'LineWidth',1);
p3 = plot(time3,x3*10^6,'LineWidth',1);
p4 = plot(timepi,xpi*10^6,'LineWidth',1);
p5 = plot(time4,x4*10^6,'LineWidth',1);
grid on
title('Displacement Function of Time for Constant Vdc')
ylabel('Displacement (um)');
xlabel('Time (s)')
legend([p5 p4 p3 p2 p1],'9.7V','Vpi=9.34V','9.3V','9V','7V','Location','best')







