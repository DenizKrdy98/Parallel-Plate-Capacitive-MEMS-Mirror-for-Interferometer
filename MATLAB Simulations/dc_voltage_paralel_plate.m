%% DC VOLTAGE Paralel Plate
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

%% Pull-in Calculation With ODE
% ODE related Values
t_f = 200e-3;
x_i = 0;
xmax = 30e-6; 
tol = 1e-24;
tol_option = odeset('AbsTol',tol);

% ODE Solution to find Vpi
Vdc = 9:0.001:12;        %***Vdc = 9:0.001:12
for i=1:length(Vdc)
    [t1,xt] = ode45(@(t,x) motion_eq2(t,x,k,Vdc(i),m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
    x_pi = xt(:,1);
    if(x_pi(end)>=g0-gmin)
        Vpi=Vdc(i)
        break;
    end
end

%% Input Voltage Sweep Configurations
Vin = 15;       %**Vin = 15, 1e4;
Vfwd = linspace(0,Vin,1e4);      % Forward Voltage
Vbwd = linspace(Vin,0,1e4);      % Backward Voltage
Vtot = [Vfwd,Vbwd];

len = length(Vfwd);             % vector length of Vforward
time = linspace(0,10,2*len);
dt = diff(time);
x = zeros(1,2*len);             % displacement 
der_x = zeros(1,2*len);         % first derivative of displacement%each element of x array
der2_x = zeros(1,2*len);        % second derivative of displacement
x_i = 0;                        % each element of x array

%%  Calculation of X (displacement)
% - Forward Voltage Calculations
check = 1;
for i=1:len
    prev_xi = x_i;
    der_x = diff(x)./dt;
    der2_x = diff(der_x)./dt(1:end-1);
    x_i = (1/k)*(Vfwd(i).^2*A*eps/(2*(g0-prev_xi)^2 - b*der_x(i) - m*der2_x(i)));
    if (Vfwd(i)>Vpi) && check==1
        Xpi = x_i*10^6
        check=0;
    end
    if x_i > g0-gmin
        x_i = g0-gmin;
    end
    x(i)=x_i;
end

% - Backward Voltage Calculations
for i=1:len
    prev_xi = x_i;
    der_x = diff(x)./dt;
    der2_x = diff(der_x)./dt(1:end-1);
    if(i<len-1)
        x_i = (1/k)*(Vbwd(i).^2*A*eps/(2*(g0-prev_xi)^2 - b*der_x(i+len) - m*der2_x(i+len)));
    end
    if x_i > g0-gmin
        x_i = g0-gmin;
    end
    x(i+len)=x_i;
end

%% Plots

% Displacement and Voltage Plot vs Time
figure
yyaxis left
plot(time,x*1e6);
ylabel('Displacement (um)');
hold on
yyaxis right
plot(time,Vtot)
ylabel('Voltage (V)');
legend('Displacement','Voltage','Location','northwest');
xlabel('Time (s)')
title('Displacement and DC Voltage vs Time in ES Actuator')
grid on

% Displacement vs Voltage Plot 
figure
plot(Vfwd,x(1:len)*1e6,'LineWidth',0.75);
hold on
plot(Vbwd,x(1+len:2*len)*1e6);
plot(Vpi,Xpi,'r*')
legend('Forward Voltage','Backward Voltage','Pull-in Voltage','Location','southeast');
ylabel('Displacement (um)');
xlabel('Voltage (V)')
title('Displacement vs Voltage in ES Actuator')
grid on




