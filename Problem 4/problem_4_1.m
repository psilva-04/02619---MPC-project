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

% % Compare h1 and h2
% 
% figure
% 
% subplot(2,2,1)
% compare(data11,G11_1,G11_2)
% title('G_{11}: F_1 \rightarrow h_1')
% 
% subplot(2,2,2)
% compare(data12,G12_1,G12_2)
% title('G_{12}: F_2 \rightarrow h_1')
% 
% subplot(2,2,3)
% compare(data21,G21_1,G21_2)
% title('G_{21}: F_1 \rightarrow h_2')
% 
% subplot(2,2,4)
% compare(data22,G22_1,G22_2)
% title('G_{22}: F_2 \rightarrow h_2')
% 
% sgtitle('First vs Second Order Models - h_1 and h_2 - 10% Step')
% 
% drawnow
% 
% axs = findall(gcf,'Type','axes');
% 
% for k = 1:length(axs)
%     axs(k).XLim = [0 1200];
%     axs(k).XTick = 0:120:1200;
%     axs(k).XTickLabel = string(0:2:20);
% end
% 
% txt = findall(gcf,'Type','text');
% 
% for k = 1:length(txt)
%     if contains(string(txt(k).String),'Time (seconds)')
%         txt(k).String = 'Time [min]';
%     end
% end
% 
% 
% % Compare h3 and h4
% 
% figure
% 
% subplot(2,2,1)
% compare(data31,G31_1,G31_2)
% title('G_{31}: F_1 \rightarrow h_3')
% 
% subplot(2,2,2)
% compare(data32,G32_1,G32_2)
% title('G_{32}: F_2 \rightarrow h_3')
% 
% subplot(2,2,3)
% compare(data41,G41_1,G41_2)
% title('G_{41}: F_1 \rightarrow h_4')
% 
% subplot(2,2,4)
% compare(data42,G42_1,G42_2)
% title('G_{42}: F_2 \rightarrow h_4')
% 
% sgtitle('First vs Second Order Models - h_3 and h_4 - 10% Step')
% 
% drawnow
% 
% axs = findall(gcf,'Type','axes');
% 
% for k = 1:length(axs)
%     axs(k).XLim = [0 1200];
%     axs(k).XTick = 0:120:1200;
%     axs(k).XTickLabel = string(0:2:20);
% end
% 
% txt = findall(gcf,'Type','text');
% 
% for k = 1:length(txt)
%     if contains(string(txt(k).String),'Time (seconds)')
%         txt(k).String = 'Time [min]';
%     end
% end

dataModels = {
    data11 G11_1 G11_2 'G_{11}: F_1 \rightarrow h_1'
    data12 G12_1 G12_2 'G_{12}: F_2 \rightarrow h_1'
    data21 G21_1 G21_2 'G_{21}: F_1 \rightarrow h_2'
    data22 G22_1 G22_2 'G_{22}: F_2 \rightarrow h_2'
    data31 G31_1 G31_2 'G_{31}: F_1 \rightarrow h_3'
    data32 G32_1 G32_2 'G_{32}: F_2 \rightarrow h_3'
    data41 G41_1 G41_2 'G_{41}: F_1 \rightarrow h_4'
    data42 G42_1 G42_2 'G_{42}: F_2 \rightarrow h_4'
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

G_delay = cell(8,1);

K_delay = zeros(8,1);
tau_delay = zeros(8,1);
theta_delay = zeros(8,1);

RMSE_1 = zeros(8,1);
RMSE_2 = zeros(8,1);
RMSE_delay = zeros(8,1);

opt = compareOptions('InitialCondition','z');

for k = 1:8

    if k == 1 || k == 5
        figure
    end

    data = dataModels{k,1};

    % First-order transfer function with delay

    G_delay{k} = tfest(data,1,0,NaN);

    K_delay(k) = dcgain(G_delay{k});

    poles_delay = pole(G_delay{k});
    tau_delay(k) = -1/real(poles_delay(1));

    theta_delay(k) = G_delay{k}.IODelay;

    % Model comparison

    [ym1,~] = compare(data,dataModels{k,2},opt);
    [ym2,~] = compare(data,dataModels{k,3},opt);
    [ymD,~] = compare(data,G_delay{k},opt);

    y = data.OutputData;
    yhat1 = ym1.OutputData;
    yhat2 = ym2.OutputData;
    yhatD = ymD.OutputData;

    idx = t_id >= t_step;

    RMSE_1(k) = sqrt(mean((y(idx)-yhat1(idx)).^2));
    RMSE_2(k) = sqrt(mean((y(idx)-yhat2(idx)).^2));
    RMSE_delay(k) = sqrt(mean((y(idx)-yhatD(idx)).^2));

    % Plots

    subplot(2,2,mod(k-1,4)+1)

    plot(t_id/60,y,'Color',[0.6 0.6 0.6])
    hold on
    plot(t_id/60,yhat1,'LineWidth',1.3)
    plot(t_id/60,yhat2,'LineWidth',1.3)
    plot(t_id/60,yhatD,'LineWidth',1.3)

    grid on
    xlim([0 20])
    xticks(0:2:20)

    title(dataModels{k,4})
    xlabel('Time [min]')
    ylabel('Normalized response')

    legend( ...
        'Simulation data', ...
        sprintf('1st order: RMSE = %.4f',RMSE_1(k)), ...
        sprintf('2nd order: RMSE = %.4f',RMSE_2(k)), ...
        sprintf('1st + delay: RMSE = %.4f',RMSE_delay(k)), ...
        'Location','best', ...
        'FontSize',7)

    if k == 4
        sgtitle('Deterministic TF Comparison - h_1 and h_2')
    elseif k == 8
        sgtitle('Deterministic TF Comparison - h_3 and h_4')
    end

end


% RMSE table

RMSETable = table( ...
    names,RMSE_1,RMSE_2,RMSE_delay, ...
    'VariableNames', ...
    {'TransferFunction','RMSE_1st','RMSE_2nd','RMSE_1st_Delay'});


% Model parameters

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

DelayTable = table( ...
    names, ...
    K_delay, ...
    tau_delay, ...
    theta_delay, ...
    'VariableNames',{'TransferFunction','Gain','Tau_s','Delay_s'});

% Save tables

filename = 'TF_Parameters.xlsx';

writetable(FirstOrderTable, filename, ...
    'Sheet', 'Deterministic 1st Order');

writetable(SecondOrderTable, filename, ...
    'Sheet', 'Deterministic 2nd Order');

writetable(DelayTable, filename, ...
    'Sheet', 'Deterministic Delay');

writetable(RMSETable, filename, ...
    'Sheet', 'Deterministic RMSE');




%% Markov parameters - deterministic model

G = {
    tf(G11_1) tf(G12_2)
    tf(G21_2) tf(G22_1)
    tf(0)     tf(G32_1)
    tf(G41_1) tf(0)
};

ny = size(G,1);
nu = size(G,2);

num = cell(ny,nu);
den = cell(ny,nu);

for i = 1:ny
    for j = 1:nu
        [num{i,j},den{i,j}] = tfdata(G{i,j},'v');
    end
end

tau = zeros(ny,nu);

Ts_markov = 4;
N = 300;

[H,tH] = mimoctf2dimpulse( ...
    num,den,tau,Ts_markov,N);


%% Plot Markov parameters

figure

for i = 1:ny

    for j = 1:nu

        subplot(ny,nu,(i-1)*nu+j)

        stem(tH/60,squeeze(H(i,j,:)),'.')
        grid on
        xlim([0 20])
        xticks(0:2:20)

        ylabel(['h_' num2str(i)])

        if i == 1
            title(['F_' num2str(j)])
        end

        if i == ny
            xlabel('Time [min]')
        end

    end

end

sgtitle('Deterministic - Discrete-Time Markov Parameters')


%% Verify Markov parameters

K_markov = sum(H,3);

K_TF = cellfun(@dcgain,G);

disp('DC gain calculated from Markov parameters:')
disp(K_markov)

disp('DC gain of identified transfer functions:')
disp(K_TF)

gainError = max(abs(K_markov(:)-K_TF(:)));

fprintf('\nMaximum DC gain difference = %.3e\n',gainError)

%% Hankel singular values - deterministic model

Nmax = 300;
tol = 1e-8;

[~,~,~,~,sH] = mimoctf2dss( ...
    num,den,tau,Ts_markov,Nmax,tol);

figure

Nsv = min(20,length(sH));

semilogy(1:Nsv,sH(1:Nsv),'o-')
grid on

if Nsv > 1
    xlim([1 Nsv])
end

xlabel('Singular value index')
ylabel('Hankel singular value')
title('Deterministic - Hankel Singular Values')

yline(tol,'--','Tolerance')


% Save Markov parameters

save('Markov_deterministic.mat', ...
    'H','tH','sH','K_markov','K_TF', ...
    'gainError','Ts_markov','N','Nmax','tol','G');


%% Physical four-state space model - deterministic

a_tank = 1.2272*ones(4,1);
A_tank = 380.1327*ones(4,1);

gamma = [0.6; 0.7];
g = 981;
rho = 1;

hs = xs./(rho*A_tank);

T = (A_tank./a_tank).*sqrt(2*hs/g);

Ac = [
    -1/T(1)   0         1/T(3)    0
     0       -1/T(2)    0         1/T(4)
     0        0        -1/T(3)    0
     0        0         0        -1/T(4)
];

Bc = rho*[
    gamma(1)    0
    0           gamma(2)
    0           1-gamma(2)
    1-gamma(1)  0
];

Ec = rho*[
    0 0
    0 0
    1 0
    0 1
];

Cy = diag(1./(rho*A_tank));
Cz = Cy(1:2,:);

Dy = zeros(4,2);
Dz = zeros(2,2);


%% Discretization - ZOH

sysc_all = ss(Ac,[Bc Ec],Cy,zeros(4,4));

sysd_all = c2d(sysc_all,Ts_markov,'zoh');

Ad = sysd_all.A;
Bd = sysd_all.B(:,1:2);
Ed = sysd_all.B(:,3:4);

sysd_y = ss(Ad,Bd,Cy,Dy,Ts_markov);
sysd = ss(Ad,Bd,Cz,Dz,Ts_markov);

Cd = Cz;
Dd = Dz;

%% Save physical state-space model

save('sysd.mat', ...
    'sysd','sysd_y','sysd_all', ...
    'Ac','Bc','Ec', ...
    'Ad','Bd','Ed', ...
    'Cy','Cz','Dy','Dz', ...
    'Cd','Dd','Ts_markov');
