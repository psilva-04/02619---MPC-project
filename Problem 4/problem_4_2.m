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

    % Save medium noise responses for TF identification

    if n == 1
    
        YnormF1_low = YnormF1;
        YnormF2_low = YnormF2;
        TOU_low = TOU;
    
    end
end

%% Transfer function identification
% 10% step and medium noise only

i = 1;

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

% First vs second order comparison - h1 and h2

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

sgtitle('Low Noise - First vs Second Order - h_1 and h_2 - 10% Step')

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


% First vs second order comparison - h3 and h4

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

sgtitle('Low Noise - First vs Second Order - h_3 and h_4 - 10% Step')

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


% Transfer function parameters

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
    'VariableNames', ...
    {'TransferFunction','Gain','Tau_s'});

SecondOrderTable = table( ...
    names, ...
    K2, ...
    tau2_1, ...
    tau2_2, ...
    'VariableNames', ...
    {'TransferFunction','Gain','Tau1_s','Tau2_s'});


fprintf('\nFIRST-ORDER MODEL PARAMETERS\n')
disp(FirstOrderTable)

fprintf('\nSECOND-ORDER MODEL PARAMETERS\n')
disp(SecondOrderTable)


% Save tables

filename = 'TF_Parameters.xlsx';

writetable(FirstOrderTable, filename, ...
    'Sheet', 'Stochastic 1st Order');

% writetable(SecondOrderTable, filename, ...
%     'Sheet', 'Stochastic 2nd Order');


%% Markov parameters - stochastic model

G_stochastic = {
    G11_1 G12_1
    G21_1 G22_1
};

num_stochastic = cell(height(G_stochastic),2);
den_stochastic = cell(height(G_stochastic),2);

for i = 1:height(G_stochastic)
    for j = 1:2
        [num_stochastic{i,j},den_stochastic{i,j}] = ...
            tfdata(G_stochastic{i,j},'v');
    end
end

tau_stochastic = zeros(height(G_stochastic),2);

Ts_markov = 4;
N = 300;

[H_stochastic,tH_stochastic] = mimoctf2dimpulse( ...
    num_stochastic,den_stochastic,tau_stochastic,Ts_markov,N);


figure

for i = 1:height(G_stochastic)

    subplot(height(G_stochastic),2,2*i-1)
    stem(tH_stochastic/60,squeeze(H_stochastic(i,1,:)),'.')
    grid on
    xlim([0 20])
    xticks(0:2:20)
    ylabel(['h_' num2str(i)])

    if i == 1
        title('F_1')
    end

    if i == height(G_stochastic)
        xlabel('Time [min]')
    end


    subplot(height(G_stochastic),2,2*i)
    stem(tH_stochastic/60,squeeze(H_stochastic(i,2,:)),'.')
    grid on
    xlim([0 20])
    xticks(0:2:20)

    if i == 1
        title('F_2')
    end

    if i == height(G_stochastic)
        xlabel('Time [min]')
    end

end

sgtitle('Low Noise - Discrete-Time Markov Parameters')


%% Discrete-time state-space realization - stochastic model

Nmax = 300;
tol = 1e-8;

[Ad_stochastic,Bd_stochastic,Cd_stochastic,Dd_stochastic,sH_stochastic] = ...
    mimoctf2dss( ...
    num_stochastic,den_stochastic,tau_stochastic, ...
    Ts_markov,Nmax,tol);

sysd_stochastic = ss(Ad_stochastic,Bd_stochastic,Cd_stochastic,Dd_stochastic,Ts_markov);
save('sysd_stochastic.mat', ...
    'sysd_stochastic', ...
    'Ad_stochastic','Bd_stochastic', ...
    'Cd_stochastic','Dd_stochastic', ...
    'Ts_markov');

% Plot Singular Values

figure

Nsv = min(20,length(sH_stochastic));

semilogy(1:Nsv,sH_stochastic(1:Nsv),'o-')
grid on
xlim([1 Nsv])

xlabel('Singular value index')
ylabel('Hankel singular value')
title('Low Noise - Hankel Singular Values')

yline(tol,'--','Tolerance')


%% Verify Markov parameters from state-space model

Hss_stochastic = zeros(height(G_stochastic),2,N+1);

Hss_stochastic(:,:,1) = Dd_stochastic;

for k = 1:N
    Hss_stochastic(:,:,k+1) = ...
        Cd_stochastic*(Ad_stochastic^(k-1))*Bd_stochastic;
end

maxError_stochastic = ...
    max(abs(H_stochastic(:)-Hss_stochastic(:)));

fprintf('\nMaximum stochastic Markov parameter realization error = %.3e\n\n', ...
    maxError_stochastic)

K_markov = sum(H_stochastic,3);

disp('DC gain calculated from Markov parameters:')
disp(K_markov)

disp('DC gain of identified transfer functions:')
disp([
    dcgain(G11_1) dcgain(G12_1)
    dcgain(G21_1) dcgain(G22_1)
])

save('Markov_Hankel_stochastic.mat', ...
    'H_stochastic','tH_stochastic','sH_stochastic', ...
    'Ts_markov','N','Nmax','tol')