%% Problem 2.2 - 2.4
% Stochastic nonlinear models and open-loop simulations.

close all
clear
clc

rng(1)

%% Parameters

a = 1.2272*ones(4,1);       % Outlet areas [cm^2]
A = 380.1327*ones(4,1);     % Tank areas [cm^2]

gamma1 = 0.6;
gamma2 = 0.7;

g = 981;                    % Gravity [cm/s^2]
rho = 1;                    % Water density [g/cm^3]

% p = [a1 a2 a3 a4 A1 A2 A3 A4 gamma1 gamma2 g rho]'
p = [a; A; gamma1; gamma2; g; rho];

%% Simulation time

dt = 0.1;
tf = 1200;

t = (0:dt:tf)';
N = length(t);

%% Piecewise-constant manipulated inputs

u = zeros(N,2);

idx = t < 400;
u(idx,:) = repmat([275 325],sum(idx),1);

idx = t >= 400 & t < 800;
u(idx,:) = repmat([300 325],sum(idx),1);

idx = t >= 800;
u(idx,:) = repmat([300 350],sum(idx),1);

%% Nominal disturbances

dMean = [50; 50];           % Mean F3 and F4 [cm^3/s]

%% Measurement noise

sigmaY = [0.10; 0.10];      % Measurement std [cm]

Rv = diag(sigmaY.^2);

%% Initial steady state

x0 = steadyStateMQT( ...
    u(1,:)', ...
    dMean, ...
    p);

%% 2.2 Piecewise-constant stochastic disturbance

% A new random disturbance value is generated every TsDist seconds
% and remains constant between updates.

TsDist = 20;                % Disturbance update period [s]

sigmaPiecewise = [3; 3];    % Disturbance std [cm^3/s]

Qd = diag(sigmaPiecewise.^2);

[xPC,yPC,zPC,dPC] = simulatePiecewise( ...
    t, ...
    u, ...
    x0, ...
    dMean, ...
    Qd, ...
    TsDist, ...
    Rv, ...
    p);

%% 2.3 Brownian disturbance

% dF = SigmaB dW

sigmaBrownian = [0.35; 0.35];

SigmaB = diag(sigmaBrownian);

[xBM,yBM,zBM,dBM] = simulateBrownian( ...
    t, ...
    u, ...
    x0, ...
    dMean, ...
    SigmaB, ...
    Rv, ...
    p);

%% 2.3 Ornstein-Uhlenbeck disturbance

% dF = Kappa*(dMean - F)dt + SigmaOU dW

kappa = [0.02; 0.02];

Kappa = diag(kappa);

sigmaOU = [1.0; 1.0];

SigmaOU = diag(sigmaOU);

[xOU,yOU,zOU,dOU] = simulateOU( ...
    t, ...
    u, ...
    x0, ...
    dMean, ...
    Kappa, ...
    SigmaOU, ...
    Rv, ...
    p);

%% Nominal steady states for each input region

uss1 = [275;325];
uss2 = [300;325];
uss3 = [300;350];

xss1 = steadyStateMQT(uss1,dMean,p);
xss2 = steadyStateMQT(uss2,dMean,p);
xss3 = steadyStateMQT(uss3,dMean,p);

zss1 = FourTankOutputs(xss1,p);
zss2 = FourTankOutputs(xss2,p);
zss3 = FourTankOutputs(xss3,p);

fprintf('\nNominal steady-state outputs\n')

fprintf('\nInput 1: F1 = %.1f, F2 = %.1f\n',uss1(1),uss1(2))
fprintf('h1 = %.3f cm\n',zss1(1))
fprintf('h2 = %.3f cm\n',zss1(2))

fprintf('\nInput 2: F1 = %.1f, F2 = %.1f\n',uss2(1),uss2(2))
fprintf('h1 = %.3f cm\n',zss2(1))
fprintf('h2 = %.3f cm\n',zss2(2))

fprintf('\nInput 3: F1 = %.1f, F2 = %.1f\n',uss3(1),uss3(2))
fprintf('h1 = %.3f cm\n',zss3(1))
fprintf('h2 = %.3f cm\n\n',zss3(2))

%% Plot manipulated inputs

figure

tiledlayout(2,1)

nexttile
stairs(t,u(:,1),'LineWidth',1.2)

ylabel('$F_1$ [cm$^3$/s]', ...
    'Interpreter','latex')

grid on

nexttile
stairs(t,u(:,2),'LineWidth',1.2)

ylabel('$F_2$ [cm$^3$/s]', ...
    'Interpreter','latex')

xlabel('Time [s]')

grid on

sgtitle('Manipulated inputs')

%% Plot F3

figure

tiledlayout(3,1)

nexttile
stairs(t,dPC(:,1),'LineWidth',1)

yline(dMean(1),'--')

ylabel('$F_3$', ...
    'Interpreter','latex')

title('Piecewise constant')

grid on

nexttile
plot(t,dBM(:,1),'LineWidth',1)

yline(dMean(1),'--')

ylabel('$F_3$', ...
    'Interpreter','latex')

title('Brownian motion')

grid on

nexttile
plot(t,dOU(:,1),'LineWidth',1)

yline(dMean(1),'--')

ylabel('$F_3$', ...
    'Interpreter','latex')

xlabel('Time [s]')

title('Ornstein-Uhlenbeck')

grid on

sgtitle('Disturbance $F_3$', ...
    'Interpreter','latex')

%% Plot F4

figure

tiledlayout(3,1)

nexttile
stairs(t,dPC(:,2),'LineWidth',1)

yline(dMean(2),'--')

ylabel('$F_4$', ...
    'Interpreter','latex')

title('Piecewise constant')

grid on

nexttile
plot(t,dBM(:,2),'LineWidth',1)

yline(dMean(2),'--')

ylabel('$F_4$', ...
    'Interpreter','latex')

title('Brownian motion')

grid on

nexttile
plot(t,dOU(:,2),'LineWidth',1)

yline(dMean(2),'--')

ylabel('$F_4$', ...
    'Interpreter','latex')

xlabel('Time [s]')

title('Ornstein-Uhlenbeck')

grid on

sgtitle('Disturbance $F_4$', ...
    'Interpreter','latex')

%% Plot h1

figure

plot(t,zPC(:,1),'LineWidth',1.2)
hold on

plot(t,zBM(:,1),'LineWidth',1.2)
plot(t,zOU(:,1),'LineWidth',1.2)

yline(zss1(1),'--')
yline(zss2(1),'--')
yline(zss3(1),'--')

xline(400,':')
xline(800,':')

xlabel('Time [s]')

ylabel('$h_1$ [cm]', ...
    'Interpreter','latex')

legend( ...
    'Piecewise constant', ...
    'Brownian', ...
    'Ornstein-Uhlenbeck', ...
    'Location','best')

title('Tank 1 output')

grid on

%% Plot h2

figure

plot(t,zPC(:,2),'LineWidth',1.2)
hold on

plot(t,zBM(:,2),'LineWidth',1.2)
plot(t,zOU(:,2),'LineWidth',1.2)

yline(zss1(2),'--')
yline(zss2(2),'--')
yline(zss3(2),'--')

xline(400,':')
xline(800,':')

xlabel('Time [s]')

ylabel('$h_2$ [cm]', ...
    'Interpreter','latex')

legend( ...
    'Piecewise constant', ...
    'Brownian', ...
    'Ornstein-Uhlenbeck', ...
    'Location','best')

title('Tank 2 output')

grid on

%% Measurement example

figure

plot(t,zBM(:,1),'LineWidth',1.2)

hold on

plot(t,yBM(:,1),'.','MarkerSize',3)

xlabel('Time [s]')

ylabel('$h_1$ [cm]', ...
    'Interpreter','latex')

legend( ...
    'True output $z_1$', ...
    'Measurement $y_1$', ...
    'Interpreter','latex', ...
    'Location','best')

title('Measurement noise')

grid on

%% Local functions

function [x,y,z,d] = simulatePiecewise( ...
    t,u,x0,dMean,Qd,TsDist,Rv,p)
% Simulates piecewise-constant Gaussian disturbances.

N = length(t);
dt = t(2)-t(1);

x = zeros(N,4);
y = zeros(N,2);
z = zeros(N,2);
d = zeros(N,2);

x(1,:) = x0';

Ld = chol(Qd,'lower');
Lv = chol(Rv,'lower');

nHold = round(TsDist/dt);

dk = dMean + Ld*randn(2,1);

for k = 1:N

    if k > 1 && mod(k-1,nHold) == 0

        dk = dMean + Ld*randn(2,1);

    end

    d(k,:) = dk';

    zk = FourTankOutputs( ...
        x(k,:)',p);

    yk = FourTankSensors( ...
        x(k,:)',p) ...
        + Lv*randn(2,1);

    z(k,:) = zk';
    y(k,:) = yk';

    if k < N

        uk = u(k,:)';

        x(k+1,:) = rk4Step( ...
            t(k), ...
            x(k,:)', ...
            uk, ...
            dk, ...
            dt, ...
            p)';

    end

end

end


function [x,y,z,d] = simulateBrownian( ...
    t,u,x0,dMean,SigmaB,Rv,p)
% Simulates Brownian disturbances.

N = length(t);
dt = t(2)-t(1);

x = zeros(N,4);
y = zeros(N,2);
z = zeros(N,2);
d = zeros(N,2);

x(1,:) = x0';
d(1,:) = dMean';

Lv = chol(Rv,'lower');

for k = 1:N

    zk = FourTankOutputs( ...
        x(k,:)',p);

    yk = FourTankSensors( ...
        x(k,:)',p) ...
        + Lv*randn(2,1);

    z(k,:) = zk';
    y(k,:) = yk';

    if k < N

        dW = sqrt(dt)*randn(2,1);

        dk = d(k,:)';

        dNext = ...
            dk ...
            + SigmaB*dW;

        d(k+1,:) = dNext';

        uk = u(k,:)';

        x(k+1,:) = rk4Step( ...
            t(k), ...
            x(k,:)', ...
            uk, ...
            dk, ...
            dt, ...
            p)';

    end

end

end


function [x,y,z,d] = simulateOU( ...
    t,u,x0,dMean,Kappa,SigmaOU,Rv,p)
% Simulates Ornstein-Uhlenbeck disturbances.

N = length(t);
dt = t(2)-t(1);

x = zeros(N,4);
y = zeros(N,2);
z = zeros(N,2);
d = zeros(N,2);

x(1,:) = x0';
d(1,:) = dMean';

Lv = chol(Rv,'lower');

for k = 1:N

    zk = FourTankOutputs( ...
        x(k,:)',p);

    yk = FourTankSensors( ...
        x(k,:)',p) ...
        + Lv*randn(2,1);

    z(k,:) = zk';
    y(k,:) = yk';

    if k < N

        dk = d(k,:)';

        dW = sqrt(dt)*randn(2,1);

        dNext = ...
            dk ...
            + Kappa*(dMean-dk)*dt ...
            + SigmaOU*dW;

        d(k+1,:) = dNext';

        uk = u(k,:)';

        x(k+1,:) = rk4Step( ...
            t(k), ...
            x(k,:)', ...
            uk, ...
            dk, ...
            dt, ...
            p)';

    end

end

end


function xNext = rk4Step( ...
    t,x,u,d,dt,p)
% One RK4 integration step.

k1 = FourTankSystem( ...
    t,x,u,d,p);

k2 = FourTankSystem( ...
    t+dt/2, ...
    x+dt*k1/2, ...
    u,d,p);

k3 = FourTankSystem( ...
    t+dt/2, ...
    x+dt*k2/2, ...
    u,d,p);

k4 = FourTankSystem( ...
    t+dt, ...
    x+dt*k3, ...
    u,d,p);

xNext = ...
    x ...
    + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

end


function xss = steadyStateMQT(u,d,p)
% Computes the steady state for constant u and d.

a = p(1:4);
A = p(5:8);

gamma1 = p(9);
gamma2 = p(10);

g = p(11);
rho = p(12);

F1 = u(1);
F2 = u(2);

F3 = d(1);
F4 = d(2);

q3 = ...
    (1-gamma2)*F2 ...
    + F3;

q4 = ...
    (1-gamma1)*F1 ...
    + F4;

q1 = ...
    gamma1*F1 ...
    + q3;

q2 = ...
    gamma2*F2 ...
    + q4;

q = [q1;q2;q3;q4];

h = ...
    (q./a).^2 ...
    /(2*g);

xss = ...
    rho*A.*h;

end