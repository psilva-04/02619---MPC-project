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

%% Simulation settings

Ts = 4;
dt = 0.1;

N = 300;
T = 0:Ts:N*Ts;

t_step = 52;

steps = [0.10 0.25 0.50];
nSteps = length(steps);

tauF = 40;
aF = 1/tauF;

processNoise = [0.02 0.05 0.10];
sensorNoise = [1 2 3];

noiseNames = {
    'Low'
    'Medium'
    'High'
};

%% Stochastic simulations

for n = 1:3

    stdF = processNoise(n)*175;
    sigmaOU = stdF*sqrt(2*aF);

    RY = sensorNoise(n)*eye(4);
    RZ = sensorNoise(n)*eye(2);

    Y_F1 = zeros(4,N+1,nSteps);
    Y_F2 = zeros(4,N+1,nSteps);

    YnormF1 = zeros(4,N+1,nSteps);
    YnormF2 = zeros(4,N+1,nSteps);

    for i = 1:nSteps

        stepSize = steps(i);

        % F1 step

        U = repmat(us,1,N);
        U(1,T(1:N) >= t_step) = us(1)*(1 + stepSize);

        rng(1)

        [tOU,XAOU,TOU,YOU,ZOU] = ...
            SimulateFourTankSDE( ...
            xs,U,Fbar,RY,RZ,p,Ts,dt,'ou',aF,sigmaOU);

        delta_u = us(1)*stepSize;

        Y_F1(:,:,i) = YOU;
        YnormF1(:,:,i) = (YOU - ys)/delta_u;

        if i == 1
            XAOU_F1 = XAOU;
            tOU_F1 = tOU;
        end

        % F2 step

        U = repmat(us,1,N);
        U(2,T(1:N) >= t_step) = us(2)*(1 + stepSize);

        rng(1)

        [tOU,XAOU,TOU,YOU,ZOU] = ...
            SimulateFourTankSDE( ...
            xs,U,Fbar,RY,RZ,p,Ts,dt,'ou',aF,sigmaOU);

        delta_u = us(2)*stepSize;

        Y_F2(:,:,i) = YOU;
        YnormF2(:,:,i) = (YOU - ys)/delta_u;

    end

    %% Raw step responses

    figure

    for j = 1:4

        allY = ys(j);

        for i = 1:nSteps
            allY = [
                allY
                squeeze(Y_F1(j,:,i))'
                squeeze(Y_F2(j,:,i))'
            ];
        end

        lower_limit = max(0,floor(min(allY)/10)*10);
        upper_limit = (floor(max(allY)/10) + 1)*10;

        subplot(4,2,2*j-1)
        hold on

        yline(ys(j),'LineWidth',1.3)

        for i = 1:nSteps
            plot( ...
                TOU/60, ...
                squeeze(Y_F1(j,:,i)), ...
                'LineWidth',1.3)
        end

        grid on
        xlim([0 20])
        xticks(0:2:20)
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
            plot( ...
                TOU/60, ...
                squeeze(Y_F2(j,:,i)), ...
                'LineWidth',1.3)
        end

        grid on
        xlim([0 20])
        xticks(0:2:20)
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

    sgtitle(sprintf( ...
        '%s Noise Step Responses: Process = %.0f%%, Sensor R = %dI', ...
        noiseNames{n}, ...
        100*processNoise(n), ...
        sensorNoise(n)))

    %% Normalized step responses

    figure

    for j = 1:4

        allY = [];

        for i = 1:nSteps
            allY = [
                allY
                squeeze(YnormF1(j,:,i))'
                squeeze(YnormF2(j,:,i))'
            ];
        end

        lower_limit = floor(min(allY)/0.1)*0.1;
        upper_limit = (floor(max(allY)/0.1) + 1)*0.1;

        subplot(4,2,2*j-1)
        hold on

        for i = 1:nSteps
            plot( ...
                TOU/60, ...
                squeeze(YnormF1(j,:,i)), ...
                'LineWidth',1.3)
        end

        grid on
        xlim([0 20])
        xticks(0:2:20)
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
            plot( ...
                TOU/60, ...
                squeeze(YnormF2(j,:,i)), ...
                'LineWidth',1.3)
        end

        grid on
        xlim([0 20])
        xticks(0:2:20)
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

    sgtitle(sprintf( ...
        'Normalized Step Responses - %s Noise', ...
        noiseNames{n}))

    % OU disturbances

    figure

    plot( ...
        tOU_F1/60, ...
        XAOU_F1(5,:), ...
        'LineWidth',1.3)

    hold on

    plot( ...
        tOU_F1/60, ...
        XAOU_F1(6,:), ...
        'LineWidth',1.3)

    yline( ...
        Fbar(1), ...
        '--', ...
        'Nominal', ...
        'HandleVisibility','off')

    grid on

    xlabel('Time [min]')
    ylabel('Flow rate [cm^3/s]')

    xlim([0 20])
    xticks(0:2:20)

    legend( ...
        'F_3','F_4', ...
        'Location','best')

    title(sprintf( ...
        'OU Disturbances - %s Noise (%.0f%%)', ...
        noiseNames{n}, ...
        100*processNoise(n)))

    % Save low noise responses for TF identification

    if n == 1
    
        YnormF1_low = YnormF1;
        YnormF2_low = YnormF2;
        TOU_low = TOU;
    
    end
end

%% Transfer function identification
% 25% step and low noise only

i = 2;

t_id = TOU_low';
u_id = double(t_id >= t_step);

% F1 -> h1, h2, h3, h4

y11 = squeeze(YnormF1_low(1,:,i))';
y21 = squeeze(YnormF1_low(2,:,i))';
y31 = squeeze(YnormF1_low(3,:,i))';
y41 = squeeze(YnormF1_low(4,:,i))';

data11 = iddata(y11,u_id,Ts);
data21 = iddata(y21,u_id,Ts);
data31 = iddata(y31,u_id,Ts);
data41 = iddata(y41,u_id,Ts);

G11_1 = tfest(data11,1,0);
G21_1 = tfest(data21,1,0);
G31_1 = tfest(data31,1,0);
G41_1 = tfest(data41,1,0);

G11_2 = tfest(data11,2,0);
G21_2 = tfest(data21,2,0);
G31_2 = tfest(data31,2,0);
G41_2 = tfest(data41,2,0);


% F2 -> h1, h2, h3, h4

y12 = squeeze(YnormF2_low(1,:,i))';
y22 = squeeze(YnormF2_low(2,:,i))';
y32 = squeeze(YnormF2_low(3,:,i))';
y42 = squeeze(YnormF2_low(4,:,i))';

data12 = iddata(y12,u_id,Ts);
data22 = iddata(y22,u_id,Ts);
data32 = iddata(y32,u_id,Ts);
data42 = iddata(y42,u_id,Ts);

G12_1 = tfest(data12,1,0);
G22_1 = tfest(data22,1,0);
G32_1 = tfest(data32,1,0);
G42_1 = tfest(data42,1,0);

G12_2 = tfest(data12,2,0);
G22_2 = tfest(data22,2,0);
G32_2 = tfest(data32,2,0);
G42_2 = tfest(data42,2,0);

% % First vs second order comparison - h1 and h2
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
% sgtitle('Low Noise - First vs Second Order - h_1 and h_2 - 25% Step')
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
% % First vs second order comparison - h3 and h4
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
% sgtitle('Low Noise - First vs Second Order - h_3 and h_4 - 25% Step')
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
        sgtitle('Stochastic TF Comparison - h_1 and h_2')
    elseif k == 8
        sgtitle('Stochastic TF Comparison - h_3 and h_4')
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
    'Sheet', 'Stochastic 1st Order');

writetable(SecondOrderTable, filename, ...
    'Sheet', 'Stochastic 2nd Order');

writetable(DelayTable, filename, ...
    'Sheet', 'Stochastic Delay');

writetable(RMSETable, filename, ...
    'Sheet', 'Stochastic RMSE');




%% Markov parameters - stochastic model

G_stochastic = {
    tf(G11_1) tf(G12_2)
    tf(G21_2) tf(G22_1)
    tf(0)     tf(G32_1)
    tf(G41_1) tf(0)
};

ny = size(G_stochastic,1);
nu = size(G_stochastic,2);

num_stochastic = cell(ny,nu);
den_stochastic = cell(ny,nu);

for i = 1:ny
    for j = 1:nu
        [num_stochastic{i,j},den_stochastic{i,j}] = ...
            tfdata(G_stochastic{i,j},'v');
    end
end

tau_stochastic = zeros(ny,nu);

Ts_markov = 4;
N = 300;

[H_stochastic,tH_stochastic] = mimoctf2dimpulse( ...
    num_stochastic,den_stochastic,tau_stochastic,Ts_markov,N);


%% Plot Markov parameters

figure

for i = 1:ny

    for j = 1:nu

        subplot(ny,nu,(i-1)*nu+j)

        stem(tH_stochastic/60,squeeze(H_stochastic(i,j,:)),'.')
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

sgtitle('Low Noise - Discrete-Time Markov Parameters')


%% Verify Markov parameters

K_markov_stochastic = sum(H_stochastic,3);

K_TF_stochastic = cellfun(@dcgain,G_stochastic);

disp('DC gain calculated from Markov parameters:')
disp(K_markov_stochastic)

disp('DC gain of identified transfer functions:')
disp(K_TF_stochastic)

gainError_stochastic = max(abs( ...
    K_markov_stochastic(:)-K_TF_stochastic(:)));

fprintf('\nMaximum stochastic DC gain difference = %.3e\n', ...
    gainError_stochastic)


%% Hankel singular values - stochastic model

Nmax = 300;
tol = 1e-8;

[~,~,~,~,sH_stochastic] = mimoctf2dss( ...
    num_stochastic,den_stochastic,tau_stochastic, ...
    Ts_markov,Nmax,tol);

figure

Nsv = min(20,length(sH_stochastic));

semilogy(1:Nsv,sH_stochastic(1:Nsv),'o-')
grid on

if Nsv > 1
    xlim([1 Nsv])
end

xlabel('Singular value index')
ylabel('Hankel singular value')
title('Low Noise - Hankel Singular Values')

yline(tol,'--','Tolerance')


%% Save Markov parameters

save('Markov_stochastic.mat', ...
    'H_stochastic','tH_stochastic','sH_stochastic', ...
    'K_markov_stochastic','K_TF_stochastic', ...
    'gainError_stochastic', ...
    'Ts_markov','N','Nmax','tol','G_stochastic');


%% Physical four-state space model - stochastic

a_tank = 1.2272*ones(4,1);
A_tank = 380.1327*ones(4,1);

gamma = [0.6; 0.7];
g = 981;
rho = 1;

hs = xs./(rho*A_tank);

T = (A_tank./a_tank).*sqrt(2*hs/g);

Ac_stochastic = [
    -1/T(1)   0         1/T(3)    0
     0       -1/T(2)    0         1/T(4)
     0        0        -1/T(3)    0
     0        0         0        -1/T(4)
];

Bc_stochastic = rho*[
    gamma(1)    0
    0           gamma(2)
    0           1-gamma(2)
    1-gamma(1)  0
];

Ec_stochastic = rho*[
    0 0
    0 0
    1 0
    0 1
];

Cy_stochastic = diag(1./(rho*A_tank));
Cz_stochastic = Cy_stochastic(1:2,:);

Dy_stochastic = zeros(4,2);
Dz_stochastic = zeros(2,2);


%% Discretization - ZOH

sysc_stochastic_all = ss( ...
    Ac_stochastic, ...
    [Bc_stochastic Ec_stochastic], ...
    Cy_stochastic,zeros(4,4));

sysd_stochastic_all = c2d( ...
    sysc_stochastic_all,Ts_markov,'zoh');

Ad_stochastic = sysd_stochastic_all.A;
Bd_stochastic = sysd_stochastic_all.B(:,1:2);
Ed_stochastic = sysd_stochastic_all.B(:,3:4);

sysd_stochastic_y = ss( ...
    Ad_stochastic,Bd_stochastic, ...
    Cy_stochastic,Dy_stochastic,Ts_markov);

sysd_stochastic = ss( ...
    Ad_stochastic,Bd_stochastic, ...
    Cz_stochastic,Dz_stochastic,Ts_markov);

Cd_stochastic = Cz_stochastic;
Dd_stochastic = Dz_stochastic;


%% Save physical state-space model

save('sysd_stochastic.mat', ...
    'sysd_stochastic', ...
    'sysd_stochastic_y', ...
    'sysd_stochastic_all', ...
    'Ac_stochastic','Bc_stochastic','Ec_stochastic', ...
    'Ad_stochastic','Bd_stochastic','Ed_stochastic', ...
    'Cy_stochastic','Cz_stochastic', ...
    'Dy_stochastic','Dz_stochastic', ...
    'Cd_stochastic','Dd_stochastic', ...
    'Ts_markov');
