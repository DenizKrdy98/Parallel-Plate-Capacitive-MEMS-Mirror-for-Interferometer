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
t_f = 1000e-3;
x_i = 0;
tol = 1e-24;
tol_option = odeset('AbsTol',tol);

%% HOW to ACHIEVE a given displacement
% Change Vdc
um=1e6;
Vac = 1;
Vdc = 1 : 0.1 : 8 ;
w=2*pi*50;
centerX = zeros(1,length(Vdc));
deltaX = zeros(1,length(Vdc));

for i=1:length(Vdc)
    
    [time,x_t] = ode45(@(t,x) motion_eq3(t,x,k,Vac,Vdc(i),w,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
    dis = x_t(:,1);
    centerX(i) = mean(dis(500:end))*um;
    Xac = max(dis(500:end)) * um - centerX(i);
    deltaX(i) = 2*Xac;
end

figure
plot(Vdc,centerX);
hold on
plot(Vdc,deltaX);
grid on
grid minor
legend('Center of Oscillation','DeltaX of Oscillation','Location','northwest');
title('For Vac=1V and fac=50Hz, Vdc vs Displacement Behavior')
ylabel('Displacement (um)')
xlabel('DC Voltage (V)')

% Change Vac
Vdc = 1;
Vac = 1 : 0.1 : 9 ;
w=2*pi*50;

centerX = zeros(1,length(Vac));
deltaX = zeros(1,length(Vac));

for i=1:length(Vac)
    
    [time,x_t] = ode45(@(t,x) motion_eq3(t,x,k,Vac(i),Vdc,w,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
    dis = x_t(:,1);
    centerX(i) = mean(dis(500:end))*um;
    Xac = max(dis(500:end)) * um - centerX(i);
    deltaX(i) = 2*Xac;
end

figure
plot(Vac,centerX);
hold on
plot(Vac,deltaX);
grid on
grid minor
legend('Center of Oscillation','DeltaX of Oscillation','Location','northwest');
title('For Vdc=1V and fac=50Hz, Vac vs Displacement Behavior')
ylabel('Displacement (um)')
xlabel('AC Voltage (V)')


%% FREQUENCY EFFECT

f1 = 50;                %frequency of voltage
w1 = 2*pi*f1;           %natural frequency of voltage
f2 = 75;               %frequency of voltage
w2 = 2*pi*f2;           %natural frequency of voltage
f3 = 100;               %frequency of voltage
w3 = 2*pi*f3;           %natural frequency of voltage
f4 = 200;                %frequency of voltage
w4 = 2*pi*f4;           %natural frequency of voltage

Vac = [1 1 1 1];
Vdc = [1 1 1 1];
[time1,x_t1] = ode45(@(t,x) motion_eq3(t,x,k,Vac(1),Vdc(1),w1,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis1 = x_t1(:,1);
[time2,x_t2] = ode45(@(t,x) motion_eq3(t,x,k,Vac(1),Vdc(2),w2,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis2 = x_t2(:,1);
[time3,x_t3] = ode45(@(t,x) motion_eq3(t,x,k,Vac(1),Vdc(3),w3,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis3 = x_t3(:,1);
[time4,x_t4] = ode45(@(t,x) motion_eq3(t,x,k,Vac(1),Vdc(4),w4,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis4 = x_t4(:,1);

V1 = Vac(1)*cos(w1*time1)+Vdc(1);
V2 = Vac(2)*cos(w1*time1)+Vdc(2);
V3 = Vac(3)*cos(w1*time1)+Vdc(3);
V4 = Vac(4)*cos(w1*time1)+Vdc(4);


um = 10^6;
x1_center_dc_eff = mean(dis1(500:end))*um
x2_center_dc_eff = mean(dis2(500:end))*um
x1_center_ac_eff = mean(dis3(500:end))*um
x2_center = mean(dis4(500:end))*um

Xac1 = max (dis1(500:end))*um - x1_center_dc_eff
Xac2 = max (dis2(500:end))*um - x2_center_dc_eff
Xac3 = max (dis3(500:end))*um - x1_center_ac_eff
Xac4 = max (dis4(500:end))*um - x2_center

figure 
plot(time1,dis1*10^6);
hold on
plot(time2,dis2*10^6);
plot(time3,dis3*10^6);
plot(time4,dis4*10^6);
grid on
grid minor
legend('fac=0.5f0','fac=0.75f0','fac=f0','fac=2f0');
xlabel('Time (s)')
ylabel('Displacement (um)')
title('Frequency Effect on Displacement Behavior (Vac=Vdc=1V)')
xlim([0.02 0.04])
 

%% Vdc and Vac Components Effects

f = 50;               %frequency of voltage
w = 2*pi*f;           %natural frequency of voltage


Vac = [1 1 2 4 ];
Vdc = [2 4 4 4 ];
[time1,x_t1] = ode45(@(t,x) motion_eq3(t,x,k,Vac(1),Vdc(1),w,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis1 = x_t1(:,1);
[time2,x_t2] = ode45(@(t,x) motion_eq3(t,x,k,Vac(2),Vdc(2),w,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis2 = x_t2(:,1);
[time3,x_t3] = ode45(@(t,x) motion_eq3(t,x,k,Vac(3),Vdc(3),w,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis3 = x_t3(:,1);
[time4,x_t4] = ode45(@(t,x) motion_eq3(t,x,k,Vac(4),Vdc(4),w,m,b,eps,A,g0) , [0,t_f], [x_i,x_i], tol_option);
dis4 = x_t4(:,1);

V1 = Vac(1)*cos(w*time1)+Vdc(1);
V2 = Vac(2)*cos(w*time1)+Vdc(2);
V3 = Vac(3)*cos(w*time1)+Vdc(3);
V4 = Vac(4)*cos(w*time1)+Vdc(4);


um = 10^6;
x1_center_dc_eff = mean(dis1(500:end))*um
x2_center_dc_eff = mean(dis2(500:end))*um
Xac1 = max (dis1(500:end))*um - x1_center_dc_eff
Xac2 = max (dis2(500:end))*um - x2_center_dc_eff

x1_center_ac_eff = mean(dis3(500:end))*um
x2_center_ac_eff = mean(dis4(500:end))*um
Xac1 = max (dis3(500:end))*um - x1_center_ac_eff
Xac2 = max (dis4(500:end))*um - x2_center

figure 
plot(time1,dis1*10^6);
hold on
plot(time2,dis2*10^6);
xlabel('Time (s)')
ylabel('Displacement (um)')
grid on
grid minor
legend('Vdc=2V','Vdc=2V');
title('Vdc Effect on Xt Behavior (Vac=1V,f0=50Hz)')
axis([0.02 0.1 0 30])

figure
plot(time3,dis3*10^6);
hold on
plot(time4,dis4*10^6);
xlabel('Time (s)')
ylabel('Displacement (um)')
grid on
grid minor
legend('Vac=2V','Vac=4V');
title('Vac Effect on Xt Behavior (Vdc=4V,f0=50Hz)')
axis([0.02 0.1 0 30])
