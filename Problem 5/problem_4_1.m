clear
close all
clc


%% Parameters

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


%% Nominal operating point

us = [
    300
    300
];

% Disturbances kept constant
Fbar = [
    175
    175
];


%% Steady state

xguess = 10000*ones(4,1);

options = optimoptions( ...
    'fsolve', ...
    'Display','none');

xs = fsolve( ...
    @(x) FourTankSystem(0,x,us,Fbar,p), ...
    xguess, ...
    options);

ys = FourTankSensors(xs,p);

zs = FourTankOutputs(xs,p);

% Simulation settings

tspan = [0 1200]; % 20 min
t_step = 50;

steps = [0.10 0.25 0.50];

nSteps = length(steps);

%% 4.1 - Deterministic nonlinear step responses

Y_F1 = cell(nSteps,1);
Y_F2 = cell(nSteps,1);

YnormF1 = cell(nSteps,1);
YnormF2 = cell(nSteps,1);

timeF1 = cell(nSteps,1);
timeF2 = cell(nSteps,1);

ys = FourTankSensors(xs,p);

for i = 1:nSteps

    stepSize = steps(i);


    % F1 step

    u_after = us;
    u_after(1) = us(1)*(1 + stepSize);

    ufun = @(t) us*(t < t_step) + u_after*(t >= t_step);

    [t,x] = ode45( ...
        @(t,x) FourTankSystem(t,x,ufun(t),Fbar,p), ...
        tspan, ...
        xs);

    y = zeros(length(t),4);

    for k = 1:length(t)
        y(k,:) = FourTankSensors(x(k,:)',p)';
    end

    delta_u = us(1)*stepSize;

    Y_F1{i} = y;
    YnormF1{i} = (y - ys')/delta_u;
    timeF1{i} = t;


    % F2 step

    u_after = us;
    u_after(2) = us(2)*(1 + stepSize);

    ufun = @(t) us*(t < t_step) + u_after*(t >= t_step);

    [t,x] = ode45( ...
        @(t,x) FourTankSystem(t,x,ufun(t),Fbar,p), ...
        tspan, ...
        xs);

    y = zeros(length(t),4);

    for k = 1:length(t)
        y(k,:) = FourTankSensors(x(k,:)',p)';
    end

    delta_u = us(2)*stepSize;

    Y_F2{i} = y;
    YnormF2{i} = (y - ys')/delta_u;
    timeF2{i} = t;

end


% Step responses

figure

for j = 1:4

    allY = ys(j);

    for i = 1:nSteps
        allY = [allY; Y_F1{i}(:,j); Y_F2{i}(:,j)];
    end

    lower_limit = max(0, floor(min(allY)/10)*10);
    upper_limit = (floor(max(allY)/10) + 1)*10;


    subplot(4,2,2*j-1)
    hold on

    yline(ys(j),'LineWidth',1.3)

    for i = 1:nSteps
        plot(timeF1{i}/60,Y_F1{i}(:,j),'LineWidth',1.3)
    end

    grid on
    ylim([lower_limit upper_limit])
    ylabel(sprintf('h_%d [cm]',j))

    if j == 1
        title('F_1')
    end

    if j == 4
        xlabel('Time [min]')
    end


    subplot(4,2,2*j)
    hold on

    yline(ys(j),'LineWidth',1.3)

    for i = 1:nSteps
        plot(timeF2{i}/60,Y_F2{i}(:,j),'LineWidth',1.3)
    end

    grid on
    ylim([lower_limit upper_limit])

    if j == 1
        title('F_2')
    end

    if j == 4
        xlabel('Time [min]')

        legend( ...
            '0% Increase', ...
            '10% Increase', ...
            '25% Increase', ...
            '50% Increase', ...
            'Location','best')
    end

end

sgtitle('Deterministic Step Responses')

% Normalized step responses

figure

for j = 1:4

    allY = [];

    for i = 1:nSteps
        allY = [allY; YnormF1{i}(:,j); YnormF2{i}(:,j)];
    end

    lower_limit = max(0, floor(min(allY)/0.1)*0.1);
    upper_limit = (floor(max(allY)/0.1) + 1)*0.1;


    subplot(4,2,2*j-1)
    hold on

    for i = 1:nSteps
        plot(timeF1{i}/60,YnormF1{i}(:,j),'LineWidth',1.3)
    end

    grid on
    ylim([lower_limit upper_limit])
    ylabel(sprintf('\\Deltah_%d / \\DeltaF_1',j))

    if j == 1
        title('F_1')
    end

    if j == 4
        xlabel('Time [min]')
    end


    subplot(4,2,2*j)
    hold on

    for i = 1:nSteps
        plot(timeF2{i}/60,YnormF2{i}(:,j),'LineWidth',1.3)
    end

    grid on
    ylim([lower_limit upper_limit])
    ylabel(sprintf('\\Deltah_%d / \\DeltaF_2',j))

    if j == 1
        title('F_2')
    end

    if j == 4
        xlabel('Time [min]')

        legend( ...
            '10% Increase', ...
            '25% Increase', ...
            '50% Increase', ...
            'Location','best')
    end

end

sgtitle('Normalized Deterministic Step Responses')

%% Transfer function identification - 10% step

Ts_id = 4;
t_start_id = 0;

i = 1;

t_id = (t_start_id:Ts_id:tspan(2))';
u_id = double(t_id >= t_step);


% F1 -> h1, h2, h3, h4

t = timeF1{i};
yNorm = YnormF1{i};

y1 = interp1(t,yNorm(:,1),t_id);
y2 = interp1(t,yNorm(:,2),t_id);
y3 = interp1(t,yNorm(:,3),t_id);
y4 = interp1(t,yNorm(:,4),t_id);

data11 = iddata(y1,u_id,Ts_id);
data21 = iddata(y2,u_id,Ts_id);
data31 = iddata(y3,u_id,Ts_id);
data41 = iddata(y4,u_id,Ts_id);

G11_1 = tfest(data11,1,0);
G21_1 = tfest(data21,1,0);
G31_1 = tfest(data31,1,0);
G41_1 = tfest(data41,1,0);

G11_2 = tfest(data11,2,0);
G21_2 = tfest(data21,2,0);
G31_2 = tfest(data31,2,0);
G41_2 = tfest(data41,2,0);


% F2 -> h1, h2, h3, h4

t = timeF2{i};
yNorm = YnormF2{i};

y1 = interp1(t,yNorm(:,1),t_id);
y2 = interp1(t,yNorm(:,2),t_id);
y3 = interp1(t,yNorm(:,3),t_id);
y4 = interp1(t,yNorm(:,4),t_id);

data12 = iddata(y1,u_id,Ts_id);
data22 = iddata(y2,u_id,Ts_id);
data32 = iddata(y3,u_id,Ts_id);
data42 = iddata(y4,u_id,Ts_id);

G12_1 = tfest(data12,1,0);
G22_1 = tfest(data22,1,0);
G32_1 = tfest(data32,1,0);
G42_1 = tfest(data42,1,0);

G12_2 = tfest(data12,2,0);
G22_2 = tfest(data22,2,0);
G32_2 = tfest(data32,2,0);
G42_2 = tfest(data42,2,0);

% Compare h1 and h2

figure

subplot(2,2,1)
compare(data11,G11_1,G11_2)
title('G_{11}: F_1 \rightarrow h_1')

subplot(2,2,2)
compare(data12,G12_1,G12_2)
title('G_{12}: F_2 \rightarrow h_1')

subplot(2,2,3)
compare(data21,G21_1,G21_2)
title('G_{21}: F_1 \rightarrow h_2')

subplot(2,2,4)
compare(data22,G22_1,G22_2)
title('G_{22}: F_2 \rightarrow h_2')

sgtitle('First vs Second Order Models - h_1 and h_2 - 10% Step')

drawnow

axs = findall(gcf,'Type','axes');

for k = 1:length(axs)
    axs(k).XLim = [0 1200];
    axs(k).XTick = 0:120:1200;
    axs(k).XTickLabel = string(0:2:20);
end

txt = findall(gcf,'Type','text');

for k = 1:length(txt)
    if contains(string(txt(k).String),'Time (seconds)')
        txt(k).String = 'Time [min]';
    end
end


% Compare h3 and h4

figure

subplot(2,2,1)
compare(data31,G31_1,G31_2)
title('G_{31}: F_1 \rightarrow h_3')

subplot(2,2,2)
compare(data32,G32_1,G32_2)
title('G_{32}: F_2 \rightarrow h_3')

subplot(2,2,3)
compare(data41,G41_1,G41_2)
title('G_{41}: F_1 \rightarrow h_4')

subplot(2,2,4)
compare(data42,G42_1,G42_2)
title('G_{42}: F_2 \rightarrow h_4')

sgtitle('First vs Second Order Models - h_3 and h_4 - 10% Step')

drawnow

axs = findall(gcf,'Type','axes');

for k = 1:length(axs)
    axs(k).XLim = [0 1200];
    axs(k).XTick = 0:120:1200;
    axs(k).XTickLabel = string(0:2:20);
end

txt = findall(gcf,'Type','text');

for k = 1:length(txt)
    if contains(string(txt(k).String),'Time (seconds)')
        txt(k).String = 'Time [min]';
    end
end


%% Model parameters

models1 = {
    G11_1
    G12_1
    G21_1
    G22_1
    G31_1
    G32_1
    G41_1
    G42_1
};

models2 = {
    G11_2
    G12_2
    G21_2
    G22_2
    G31_2
    G32_2
    G41_2
    G42_2
};

names = {
    'G11'
    'G12'
    'G21'
    'G22'
    'G31'
    'G32'
    'G41'
    'G42'
};

K1 = zeros(8,1);
tau1 = zeros(8,1);

K2 = zeros(8,1);
tau2_1 = zeros(8,1);
tau2_2 = zeros(8,1);

for k = 1:8

    K1(k) = dcgain(models1{k});

    p1 = pole(models1{k});
    tau1(k) = -1/real(p1(1));

    K2(k) = dcgain(models2{k});

    p2 = pole(models2{k});

    tau2_1(k) = -1/real(p2(1));
    tau2_2(k) = -1/real(p2(2));

end


FirstOrderTable = table( ...
    names, ...
    K1, ...
    tau1, ...
    'VariableNames',{'TransferFunction','Gain','Tau_s'});

SecondOrderTable = table( ...
    names, ...
    K2, ...
    tau2_1, ...
    tau2_2, ...
    'VariableNames',{'TransferFunction','Gain','Tau1_s','Tau2_s'});


fprintf('\nFIRST-ORDER MODEL PARAMETERS\n')
disp(FirstOrderTable)

fprintf('\nSECOND-ORDER MODEL PARAMETERS\n')
disp(SecondOrderTable)

% Save tables

filename = 'TF_Parameters.xlsx';

writetable(FirstOrderTable, filename, ...
    'Sheet', 'Deterministic 1st Order');

% writetable(SecondOrderTable, filename, ...
%     'Sheet', 'Deterministic 2nd Order');


%% Markov parameters

G = {
    G11_1 G12_1
    G21_1 G22_1
};

num = cell(height(G),2);
den = cell(height(G),2);

for i = 1:height(G)
    for j = 1:2
        [num{i,j},den{i,j}] = tfdata(G{i,j},'v');
    end
end

tau = zeros(height(G),2);

Ts_markov = 4;
N = 300;

[H,tH] = mimoctf2dimpulse( ...
    num,den,tau,Ts_markov,N);


figure

for i = 1:height(G)

    subplot(height(G),2,2*i-1)
    stem(tH/60,squeeze(H(i,1,:)),'.')
    grid on
    xlim([0 20])
    xticks(0:2:20)
    ylabel(['h_' num2str(i)])

    if i == 1
        title('F_1')
    end

    if i == height(G)
        xlabel('Time [min]')
    end


    subplot(height(G),2,2*i)
    stem(tH/60,squeeze(H(i,2,:)),'.')
    grid on
    xlim([0 20])
    xticks(0:2:20)

    if i == 1
        title('F_2')
    end

    if i == height(G)
        xlabel('Time [min]')
    end

end

sgtitle('Discrete-Time Markov Parameters')

%% Discrete-time state-space realization

Nmax = 300;
tol = 1e-8;

[Ad,Bd,Cd,Dd,sH] = mimoctf2dss( ...
    num,den,tau,Ts_markov,Nmax,tol);

sysd = ss(Ad,Bd,Cd,Dd,Ts_markov);
save('sysd.mat','sysd','Ad','Bd','Cd','Dd','Ts_markov');

% Plot Singular Values

figure

Nsv = min(20,length(sH));

semilogy(1:Nsv,sH(1:Nsv),'o-')
grid on
xlim([1 Nsv])

xlabel('Singular value index')
ylabel('Hankel singular value')
title('Hankel Singular Values')

yline(tol,'--','Tolerance')

%% Verify Markov parameters from state-space model

Hss = zeros(height(G),2,N+1);

Hss(:,:,1) = Dd;

for k = 1:N
    Hss(:,:,k+1) = Cd*(Ad^(k-1))*Bd;
end

maxError = max(abs(H(:)-Hss(:)));

fprintf('\nMaximum Markov parameter realization error = %.3e\n\n', ...
    maxError)

K_markov = sum(H,3);

disp('DC gain calculated from Markov parameters:')
disp(K_markov)

disp('DC gain of identified transfer functions:')
disp([
    dcgain(G11_1) dcgain(G12_1)
    dcgain(G21_1) dcgain(G22_1)
])

save('Markov_Hankel.mat', ...
    'H','tH','sH', ...
    'Ts_markov','N','Nmax','tol')