clear
close all
clc

rng(1)


%% Parameters

a = 1.2272*ones(4,1);
A = 380.1327*ones(4,1);

gamma = [0.6; 0.7];
g = 981;
rho = 1;

p = [a; A; gamma; g; rho];


%% Nominal operating point

us = [300; 300];
ds = [175; 175];


%% Steady state

xguess = 10000*ones(4,1);

options = optimoptions('fsolve','Display','none');

xs = fsolve( ...
    @(x) FourTankSystem(0,x,us,ds,p), ...
    xguess, ...
    options);

ys = FourTankSensors(xs,p);
zs = FourTankOutputs(xs,p);

disp('Steady-state masses [g]:')
disp(xs)

disp('Steady-state heights [cm]:')
disp(ys)


%% Noise covariance matrices

Qd = 306.25*eye(2);
R  = 3*eye(4);


%% Simulation settings

Ts = 5;
N = 120;

U = repmat(us,1,N);


%% Simulation

[T,X,Y,Z,D] = SimulateFourTankPWC( ...
    xs,U,ds,Qd,R,p,Ts);


%% Plot disturbances F3 and F4

figure

stairs(T(1:end-1),D(1,:),'LineWidth',1.5)
hold on
stairs(T(1:end-1),D(2,:),'LineWidth',1.5)

yline(ds(1),'--')

grid on

xlabel('Time [s]')
ylabel('Flow [cm^3/s]')
title('Unmeasured disturbances')

legend('F_3','F_4','Mean','Location','best')


%% Plot measured outputs h1 and h2

figure

plot(T,Y(1,:),'LineWidth',1.2)
hold on
plot(T,Y(2,:),'LineWidth',1.2)

grid on

xlabel('Time [s]')
ylabel('Level [cm]')
title('Measured outputs')

legend('Measured h_1','Measured h_2','Location','best')