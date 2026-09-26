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


% Noise parameters


Qd = 306.25*eye(2);

R = 3*eye(4);


% Simulation settings

Ts = 4;

dt = 0.1;

N = 80;

T = 0:Ts:N*Ts;


% Piecewise-constant manipulated inputs
%
% Four plateaus, each lasting:
%
% 20 samples * 4 s = 80 s so the system has time to respond


U = zeros(2,N);


% Plateau 1: 0 - 80 s

U(:,1:20) = repmat( ...
    [300;300], ...
    1,20);


% Plateau 2: 80 - 160 s

U(:,21:40) = repmat( ...
    [340;260], ...
    1,20);


% Plateau 3: 160 - 240 s

U(:,41:60) = repmat( ...
    [270;340], ...
    1,20);


% Plateau 4: 240 - 320 s

U(:,61:80) = repmat( ...
    [320;290], ...
    1,20);



% 2.2 - Piecewise-constant stochastic disturbances

rng(1)

[TPWC,XPWC,YPWC,ZPWC,DPWC] = ...
    SimulateFourTankPWC( ...
    xs, ...
    U, ...
    Fbar, ...
    Qd, ...
    R, ...
    p, ...
    Ts);



% 2.3 - Brownian disturbance model


sigmaB = 17.5/sqrt(Ts);

rng(1)

[tB,XAB,TB,YB,ZB] = ...
    SimulateFourTankSDE( ...
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

DB = XAB(5:6,:);


% 2.3 - OU disturbance model


tauF = 60;

aF = 1/tauF;

sigmaOU = 17.5*sqrt(2*aF);

rng(1)

[tOU,XAOU,TOU,YOU,ZOU] = ...
    SimulateFourTankSDE( ...
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

DOU = XAOU(5:6,:);



% Plot manipulated inputs F1 and F2


figure

subplot(2,1,1)

stairs(T(1:end-1),U(1,:), ...
    'LineWidth',1.5)

xlabel('Time [s]')
ylabel('F_1 [cm^3/s]')
title('Piecewise-constant input F_1')

ylim([240 360])

grid on


subplot(2,1,2)

stairs(T(1:end-1),U(2,:), ...
    'LineWidth',1.5)

xlabel('Time [s]')
ylabel('F_2 [cm^3/s]')
title('Piecewise-constant input F_2')

ylim([240 360])

grid on



% Compare true controlled output h1


figure

plot(TPWC,ZPWC(1,:), ...
    'LineWidth',1.3)

hold on

plot(TB,ZB(1,:), ...
    'LineWidth',1.3)

plot(TOU,ZOU(1,:), ...
    'LineWidth',1.3)

yline(zs(1),'--')

xlabel('Time [s]')
ylabel('h_1 [cm]')

title('Controlled output h_1')

legend( ...
    'PWC disturbance', ...
    'Brownian disturbance', ...
    'OU disturbance', ...
    'Initial steady state', ...
    'Location','best')

grid on


% Compare true controlled output h2

figure

plot(TPWC,ZPWC(2,:), ...
    'LineWidth',1.3)

hold on

plot(TB,ZB(2,:), ...
    'LineWidth',1.3)

plot(TOU,ZOU(2,:), ...
    'LineWidth',1.3)

yline(zs(2),'--')

xlabel('Time [s]')
ylabel('h_2 [cm]')

title('Controlled output h_2')

legend( ...
    'PWC disturbance', ...
    'Brownian disturbance', ...
    'OU disturbance', ...
    'Initial steady state', ...
    'Location','best')

grid on


% PWC disturbances F3 and F4


figure

subplot(2,1,1)

stairs(TPWC(1:end-1),DPWC(1,:), ...
    'LineWidth',1.2)

hold on

yline(Fbar(1),'--')

xlabel('Time [s]')
ylabel('F_3 [cm^3/s]')
title('PWC disturbance F_3')

grid on


subplot(2,1,2)

stairs(TPWC(1:end-1),DPWC(2,:), ...
    'LineWidth',1.2)

hold on

yline(Fbar(2),'--')

xlabel('Time [s]')
ylabel('F_4 [cm^3/s]')
title('PWC disturbance F_4')

grid on


% Brownian disturbances F3 and F4

figure

subplot(2,1,1)

plot(tB,DB(1,:), ...
    'LineWidth',1.2)

hold on

yline(Fbar(1),'--')

xlabel('Time [s]')
ylabel('F_3 [cm^3/s]')
title('Brownian disturbance F_3')

grid on


subplot(2,1,2)

plot(tB,DB(2,:), ...
    'LineWidth',1.2)

hold on

yline(Fbar(2),'--')

xlabel('Time [s]')
ylabel('F_4 [cm^3/s]')
title('Brownian disturbance F_4')

grid on



% OU disturbances F3 and F4

figure

subplot(2,1,1)

plot(tOU,DOU(1,:), ...
    'LineWidth',1.2)

hold on

yline(Fbar(1),'--')

xlabel('Time [s]')
ylabel('F_3 [cm^3/s]')
title('OU disturbance F_3')

grid on


subplot(2,1,2)

plot(tOU,DOU(2,:), ...
    'LineWidth',1.2)

hold on

yline(Fbar(2),'--')

xlabel('Time [s]')
ylabel('F_4 [cm^3/s]')
title('OU disturbance F_4')

grid on