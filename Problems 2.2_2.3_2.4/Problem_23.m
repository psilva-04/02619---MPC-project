clear
close all
clc


% Parameters


a = 1.2272*ones(4,1);

A = 380.1327*ones(4,1);

gamma = [
    0.6
    0.7
];

g = 981;

rho = 1;

p = [
    a
    A
    gamma
    g
    rho
];


% Nominal operating point

us = [
    300
    300
];

Fbar = [
    175
    175
];


% Steady state

xguess = 10000*ones(4,1);

options = optimoptions( ...
    'fsolve', ...
    'Display','none');

xs = fsolve( ...
    @(x) FourTankSystem(0,x,us,Fbar,p), ...
    xguess, ...
    options);


zs = FourTankOutputs(xs,p);


disp('Steady-state masses [g]:')
disp(xs)

disp('Steady-state heights [cm]:')
disp(FourTankSensors(xs,p))


% Measurement noise

R = 3*eye(4);



% Simulation settings

Ts = 5;       % Sampling time [s]

dt = 0.1;     % Euler-Maruyama integration step [s]

N = 150;       % Number of sampling intervals


% Piecewise-constant manipulated inputs

U = repmat(us,1,N);



% Brownian disturbance model
%
% dF = sigmaB dW
%
% i chose sigmaB so that the standard deviation of the  Brownian increment over Ts seconds is 17.5.

sigmaB = 17.5/sqrt(Ts);


rng(1)

[tB,XAB,T,YB,ZB] = SimulateFourTankSDE( ...
    xs, ...
    U, ...
    Fbar, ...
    R, ...
    p, ...
    Ts, ...
    dt, ...
    'brownian', ...
    0, ...
    sigmaB);


%% Separate tank states and disturbances

XB = XAB(1:4,:);

DB = XAB(5:6,:);



% Ornstein-Uhlenbeck disturbance model (which we are using from now on)
%
% dF = aF*(Fbar-F)dt + sigmaOU dW


tauF = 60;

aF = 1/tauF;


% Stationary standard deviation = 17.5

sigmaOU = 17.5*sqrt(2*aF);


rng(1)

[tOU,XAOU,~,YOU,ZOU] = SimulateFourTankSDE( ...
    xs, ...
    U, ...
    Fbar, ...
    R, ...
    p, ...
    Ts, ...
    dt, ...
    'ou', ...
    aF, ...
    sigmaOU);


%% Separate tank states and disturbances

XOU = XAOU(1:4,:);

DOU = XAOU(5:6,:);


% Brownian disturbances

figure

subplot(2,1,1)

plot(tB,DB(1,:),'LineWidth',1.2)
hold on

yline(Fbar(1),'--')

xlabel('Time [s]')
ylabel('F_3 [cm^3/s]')
title('Brownian disturbance F_3')

grid on


subplot(2,1,2)

plot(tB,DB(2,:),'LineWidth',1.2)
hold on

yline(Fbar(2),'--')

xlabel('Time [s]')
ylabel('F_4 [cm^3/s]')
title('Brownian disturbance F_4')

grid on


% OU disturbances


figure

subplot(2,1,1)

plot(tOU,DOU(1,:),'LineWidth',1.2)
hold on

yline(Fbar(1),'--')

xlabel('Time [s]')
ylabel('F_3 [cm^3/s]')
title('OU disturbance F_3')

grid on


subplot(2,1,2)

plot(tOU,DOU(2,:),'LineWidth',1.2)
hold on

yline(Fbar(2),'--')

xlabel('Time [s]')
ylabel('F_4 [cm^3/s]')
title('OU disturbance F_4')

grid on



% Measured outputs - Brownian


figure

subplot(2,1,1)

plot(T,YB(1,:),'LineWidth',1.2)
hold on

yline(zs(1),'--')

xlabel('Time [s]')
ylabel('h_1 [cm]')
title('Measured h_1 - Brownian')

grid on


subplot(2,1,2)

plot(T,YB(2,:),'LineWidth',1.2)
hold on

yline(zs(2),'--')

xlabel('Time [s]')
ylabel('h_2 [cm]')
title('Measured h_2 - Brownian')

grid on


% Measured outputs - OU


figure

subplot(2,1,1)

plot(T,YOU(1,:),'LineWidth',1.2)
hold on

yline(zs(1),'--')

xlabel('Time [s]')
ylabel('h_1 [cm]')
title('Measured h_1 - OU')

grid on


subplot(2,1,2)

plot(T,YOU(2,:),'LineWidth',1.2)
hold on

yline(zs(2),'--')

xlabel('Time [s]')
ylabel('h_2 [cm]')
title('Measured h_2 - OU')

grid on


% Comparison of true controlled outputs


figure

subplot(2,1,1)

plot(T,ZB(1,:),'LineWidth',1.2)
hold on

plot(T,ZOU(1,:),'LineWidth',1.2)

yline(zs(1),'--')

xlabel('Time [s]')
ylabel('h_1 [cm]')
title('Controlled output h_1')

legend('Brownian','OU','Nominal')

grid on


subplot(2,1,2)

plot(T,ZB(2,:),'LineWidth',1.2)
hold on

plot(T,ZOU(2,:),'LineWidth',1.2)

yline(zs(2),'--')

xlabel('Time [s]')
ylabel('h_2 [cm]')
title('Controlled output h_2')

legend('Brownian','OU','Nominal')

grid on



% FIGURE 2

figure(6)
clf
set(gcf,'Color','w','Position',[100 100 900 700])


subplot(2,1,1)

plot(tB,DB(1,:),'LineWidth',1.8)
hold on

plot(tB,DB(2,:),'LineWidth',1.8)

yline(Fbar(1),'--k','LineWidth',1.3)

xlabel('Time [s]','FontSize',14)
ylabel('DVs [cm^3/s]','FontSize',14)
title('Brownian disturbances','FontSize',16)

legend('$F_3$','$F_4$','Nominal', ...
    'Interpreter','latex', ...
    'FontSize',13, ...
    'Location','best')

grid on
box on

set(gca,'FontSize',13,'LineWidth',1)



subplot(2,1,2)

plot(T,YB(1,:),'LineWidth',1.8)
hold on

plot(T,YB(2,:),'LineWidth',1.8)

yline(zs(1),'--','LineWidth',1.3)
yline(zs(2),'--','LineWidth',1.3)

xlabel('Time [s]','FontSize',14)
ylabel('Measured Heights [cm]','FontSize',14)
title('Measured outputs - Brownian','FontSize',16)

legend('$h_1$','$h_2$', ...
    '$h_{1,s}$','$h_{2,s}$', ...
    'Interpreter','latex', ...
    'FontSize',13, ...
    'Location','best')

grid on
box on

set(gca,'FontSize',13,'LineWidth',1)


%% Brownian F3/F4 statistics

F3_B_min = min(DB(1,:));
F3_B_max = max(DB(1,:));
F3_B_range = F3_B_max - F3_B_min;

F4_B_min = min(DB(2,:));
F4_B_max = max(DB(2,:));
F4_B_range = F4_B_max - F4_B_min;

fprintf('\n========================================\n')
fprintf('BROWNIAN DISTURBANCES\n')
fprintf('========================================\n')

fprintf('F3:\n')
fprintf('  Minimum = %.3f cm^3/s\n',F3_B_min)
fprintf('  Maximum = %.3f cm^3/s\n',F3_B_max)
fprintf('  Range   = %.3f cm^3/s\n',F3_B_range)

fprintf('\nF4:\n')
fprintf('  Minimum = %.3f cm^3/s\n',F4_B_min)
fprintf('  Maximum = %.3f cm^3/s\n',F4_B_max)
fprintf('  Range   = %.3f cm^3/s\n',F4_B_range)

fprintf('========================================\n')



% FIGURE 3
% Brownian vs OU: F3 + measured h1


figure(76)
clf
set(gcf,'Color','w','Position',[150 150 900 700])


subplot(2,1,1)

plot(tB,DB(1,:),'LineWidth',1.8)
hold on

plot(tOU,DOU(1,:),'LineWidth',1.8)

yline(Fbar(1),'--k','LineWidth',1.3)

xlabel('Time [s]','FontSize',14)

ylabel('$F_3$ [cm$^3$/s]', ...
    'Interpreter','latex', ...
    'FontSize',14)

title('$F_3$: Brownian vs. Ornstein--Uhlenbeck', ...
    'Interpreter','latex', ...
    'FontSize',16)

legend('Brownian','OU','Nominal', ...
    'FontSize',13, ...
    'Location','best')

grid on
box on

set(gca,'FontSize',13,'LineWidth',1)


subplot(2,1,2)

plot(T,YB(1,:),'LineWidth',1.8)
hold on

plot(T,YOU(1,:),'LineWidth',1.8)

yline(zs(1),'--k','LineWidth',1.3)

xlabel('Time [s]','FontSize',14)

ylabel('$h_1$ [cm]', ...
    'Interpreter','latex', ...
    'FontSize',14)

title('Measured $h_1$: Brownian vs. Ornstein--Uhlenbeck', ...
    'Interpreter','latex', ...
    'FontSize',16)

legend('Brownian','OU','Nominal', ...
    'FontSize',13, ...
    'Location','best')

grid on
box on

set(gca,'FontSize',13,'LineWidth',1)


%% F3 Brownian vs OU statistics

% Recalculate here so this section also works independently

F3_B_min = min(DB(1,:));
F3_B_max = max(DB(1,:));
F3_B_range = F3_B_max - F3_B_min;

F3_OU_min = min(DOU(1,:));
F3_OU_max = max(DOU(1,:));
F3_OU_range = F3_OU_max - F3_OU_min;

fprintf('\n========================================\n')
fprintf('F3 COMPARISON: BROWNIAN vs OU\n')
fprintf('========================================\n')

fprintf('Brownian F3:\n')
fprintf('  Minimum = %.3f cm^3/s\n',F3_B_min)
fprintf('  Maximum = %.3f cm^3/s\n',F3_B_max)
fprintf('  Range   = %.3f cm^3/s\n',F3_B_range)

fprintf('\nOU F3:\n')
fprintf('  Minimum = %.3f cm^3/s\n',F3_OU_min)
fprintf('  Maximum = %.3f cm^3/s\n',F3_OU_max)
fprintf('  Range   = %.3f cm^3/s\n',F3_OU_range)

fprintf('========================================\n')


exportgraphics(figure(6), ...
    'figure6.png', ...
    'Resolution',600)

exportgraphics(figure(76), ...
    'figure76.png', ...
    'Resolution',600)