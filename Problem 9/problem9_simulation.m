clearvars; clc; close all;

% Parameters
a = 1.2272 * ones(4,1);
A = 380.1327 * ones(4,1);
gamma = [0.6; 0.7];
g = 981;
rho = 1;
p = [a; A; gamma; g; rho];

us = [300; 300];
Fbar = [175; 175];

% Nominal equilibrium masses
q3s = (1-gamma(2))*us(2) + Fbar(1);
q4s = (1-gamma(1))*us(1) + Fbar(2);

qs = [
    gamma(1)*us(1)+q3s
    gamma(2)*us(2)+q4s
    q3s
    q4s
];

hs = (qs ./ a).^2 / (2*g);
xs = rho * A .* hs;

% Load physical mass-deviation model
S = load('sysd_stochastic.mat', ...
    'Ad_stochastic', 'Bd_stochastic', 'Ed_stochastic', ...
    'Cz_stochastic', 'Ts_markov');

Ad = S.Ad_stochastic;
Bd = S.Bd_stochastic;
Ed = S.Ed_stochastic;
Cz = S.Cz_stochastic;
Ts = S.Ts_markov;

% Integration step
dt = 0.1;
nSub = round(Ts/dt);

assert(abs(nSub*dt-Ts) < 1e-10, ...
    'Ts must be an integer multiple of dt.');

nx = size(Ad,1);
nu = size(Bd,2);
nz = size(Cz,1);
nw = size(Ed,2);

Hmass = diag(1 ./ (rho*A));
Hz = Hmass(1:2,:);

assert(nx == 4 && norm(Cz-Hz,'fro') < 1e-10, ...
    'The controller must use physical mass-deviation coordinates.');

assert(norm(FourTankSystem(0,xs,us,Fbar,p),inf) < 1e-6, ...
    'The nominal masses do not satisfy the steady-state balances.');

% Controller design
Tpred = 200;
N = round(Tpred/Ts);

Wz = diag([10, 10]);
Wu = diag([0.1, 0.1]);
Wdu = eye(2);

umin = -300 * ones(nu,1);
umax = 200 * ones(nu,1);
dumin = -10 * ones(nu,1);
dumax = 10 * ones(nu,1);

mpc = mpc_design(Ad,Bd,Ed,Cz,N,Wz,Wu,Wdu, ...
    umin,umax,dumin,dumax);

% OU settings
tauF = 40;
aF = 1/tauF;
stdF = 0.05*Fbar(1);
sigmaOU = stdF * sqrt(2*aF);

alphaFine = 1-aF*dt;
alphaOU = alphaFine^nSub;

% Conditional OU prediction consistent with Euler-Maruyama
Sdist = zeros(nx,nw);
mpc.Phi_w = zeros(N*nz,nw);

for i = 1:N
    Sdist = Ad*Sdist + Ed*alphaOU^(i-1);
    rows = (i-1)*nz + (1:nz);
    mpc.Phi_w(rows,:) = Cz*Sdist;
end

% One optimisation and its next-output prediction
r = [10; -15];
R_k = repmat(r,N,1);
U_bar = zeros(N*nu,1);

x0 = zeros(nx,1);
w0 = zeros(nw,1);
u0 = zeros(nu,1);
xinit = zeros(N*nu,1);

[uk_single,Uk_single,Zk_single,info_single] = mpc_compute( ...
    mpc,x0,w0,u0,R_k,U_bar,xinit);

U_pred = reshape(Uk_single,nu,N);
du_pred = diff([u0,U_pred],1,2);
Z_pred = reshape(Zk_single,nz,N);
z_next = Z_pred(:,1);

disp('Input to apply now, deviation [cm^3/s]:');
disp(uk_single);

disp('Next predicted level deviation [cm]:');
disp(z_next);

z_check = Cz*(Ad*x0 + Bd*uk_single + Ed*w0);

fprintf('One-step consistency error: %.3e\n', ...
    norm(z_next-z_check));

figure('Name','Problem 9: single prediction');

t_pred_u = (0:N)*Ts;

subplot(3,1,1);
plot([0,(1:N)*Ts],[Cz*x0,Z_pred].','LineWidth',1.2);
hold on;
yline(r(1),'--','HandleVisibility','off');
yline(r(2),'--','HandleVisibility','off');
grid on;
xlim([0 Tpred]);
ylabel('Level deviation [cm]');
legend('Tank 1','Tank 2','Location','best');

subplot(3,1,2);
stairs(t_pred_u,[U_pred,U_pred(:,end)].','LineWidth',1.2);
hold on;
yline(umin(1),':','HandleVisibility','off');
yline(umax(1),':','HandleVisibility','off');
grid on;
xlim([0 Tpred]);
ylabel('Input deviation [cm^3/s]');
legend('Pump 1','Pump 2','Location','best');

subplot(3,1,3);
stairs(t_pred_u,[du_pred,du_pred(:,end)].','LineWidth',1.2);
hold on;
yline(dumin(1),':','HandleVisibility','off');
yline(dumax(1),':','HandleVisibility','off');
grid on;
xlim([0 Tpred]);
xlabel('Prediction time [s]');
ylabel('Input change per sample');
legend('Pump 1','Pump 2','Location','best');

% PROBLEM 11

% Shared Wiener increments and Euler-Maruyama OU realisation
Tsim = 1200;
Nsim = round(Tsim/Ts);

rng(1);

dW_fine = sqrt(dt)*randn(nw,Nsim*nSub);

w_fine = zeros(nw,Nsim*nSub+1);
w_fine(:,1) = w0;

for j = 1:Nsim*nSub
    w_fine(:,j+1) = w_fine(:,j) ...
        - aF*w_fine(:,j)*dt ...
        + sigmaOU*dW_fine(:,j);
end

w_sampled = w_fine(:,1:nSub:end);

% Independent closed-loop tests
x_linear = zeros(nx,Nsim+1);
m_nonlinear = zeros(4,Nsim+1);
m_nonlinear(:,1) = xs;

u_linear = zeros(nu,Nsim);
u_nonlinear = zeros(nu,Nsim);

z_next_linear = zeros(nz,Nsim);
z_next_nonlinear = zeros(nz,Nsim);

exitflag_linear = zeros(1,Nsim);
exitflag_nonlinear = zeros(1,Nsim);

u_prev_linear = u0;
u_prev_nonlinear = u0;

xinit_linear = zeros(N*nu,1);
xinit_nonlinear = zeros(N*nu,1);

for k = 1:Nsim
    w_k = w_sampled(:,k);

    [uL,UL,ZL,infoL] = mpc_compute( ...
        mpc,x_linear(:,k),w_k,u_prev_linear, ...
        R_k,U_bar,xinit_linear);

    x_nonlinear = m_nonlinear(:,k)-xs;

    [uN,UN,ZN,infoN] = mpc_compute( ...
        mpc,x_nonlinear,w_k,u_prev_nonlinear, ...
        R_k,U_bar,xinit_nonlinear);

    u_linear(:,k) = uL;
    u_nonlinear(:,k) = uN;

    z_next_linear(:,k) = ZL(1:nz);
    z_next_nonlinear(:,k) = ZN(1:nz);

    exitflag_linear(k) = infoL.exitflag;
    exitflag_nonlinear(k) = infoN.exitflag;

    % Discrete linear plant
    x_linear(:,k+1) = Ad*x_linear(:,k) + Bd*uL + Ed*w_k;

    % Nonlinear plant using Euler-Maruyama
    m = m_nonlinear(:,k);
    u_actual = us+uN;
    j0 = (k-1)*nSub+1;

    for j = j0:j0+nSub-1
        xa = [m; Fbar+w_fine(:,j)];

        [f,sigma] = FourTankSDE((j-1)*dt,xa,u_actual, ...
            Fbar,aF,sigmaOU,p,'ou');

        xa_next = xa + f*dt + sigma*dW_fine(:,j);
        m = xa_next(1:4);

        if ~isreal(m) || any(~isfinite(m)) || any(m < 0)
            error('Invalid nonlinear tank mass at t = %.2f s.',j*dt);
        end
    end

    m_nonlinear(:,k+1) = m;

    u_prev_linear = uL;
    u_prev_nonlinear = uN;

    xinit_linear = [UL(nu+1:end); UL(end-nu+1:end)];
    xinit_nonlinear = [UN(nu+1:end); UN(end-nu+1:end)];
end

z_linear = Cz*x_linear;
z_nonlinear = Hz*(m_nonlinear-xs);

du_linear = diff([u0,u_linear],1,2);
du_nonlinear = diff([u0,u_nonlinear],1,2);

t_x = (0:Nsim)*Ts;
t_u = (0:Nsim-1)*Ts;

figure('Name','Problem 9: closed-loop outputs');

for i = 1:nz
    subplot(nz,1,i);
    plot(t_x,z_linear(i,:),t_x,z_nonlinear(i,:), ...
        'LineWidth',1.2);
    hold on;
    yline(r(i),'--','HandleVisibility','off');
    grid on;
    ylabel(sprintf('Delta h_%d [cm]',i));
    legend('Linear plant','Nonlinear plant','Location','best');
end

xlabel('Time [s]');

figure('Name','Problem 9: closed-loop inputs');

for i = 1:nu
    subplot(nu,2,2*i-1);
    stairs(t_u,u_linear(i,:),'LineWidth',1.2);
    hold on;
    stairs(t_u,u_nonlinear(i,:),'LineWidth',1.2);
    yline(umin(i),':','HandleVisibility','off');
    yline(umax(i),':','HandleVisibility','off');
    grid on;
    ylabel(sprintf('Delta F_%d [cm^3/s]',i));
    xlabel('Time [s]');
    legend('Linear plant','Nonlinear plant','Location','best');

    subplot(nu,2,2*i);
    stairs(t_u,du_linear(i,:),'LineWidth',1.2);
    hold on;
    stairs(t_u,du_nonlinear(i,:),'LineWidth',1.2);
    yline(dumin(i),':','HandleVisibility','off');
    yline(dumax(i),':','HandleVisibility','off');
    grid on;
    ylabel(sprintf('Input %d change per sample',i));
    xlabel('Time [s]');
    legend('Linear plant','Nonlinear plant','Location','best');
end

figure('Name','Problem 9: OU disturbances');
plot((0:Nsim*nSub)*dt,(Fbar+w_fine).','LineWidth',1.2);
hold on;
yline(Fbar(1),'--','HandleVisibility','off');
grid on;
xlabel('Time [s]');
ylabel('Disturbance flow [cm^3/s]');
legend('F_3','F_4','Location','best');