% Problems 12.1 and 12.2: nonlinear plant simulation and CD-EKF.
% Run this script to reproduce the results and save the figures.
% The operating-point calculation uses fsolve (Optimization Toolbox).

clear
close all
clc

%% Parameters and nominal operating point

a = 1.2272*ones(4,1);
A_tank = 380.1327*ones(4,1);
gamma = [0.6;0.7];
g = 981;
rho = 1;
p = [a;A_tank;gamma;g;rho];

us = [300;300];          % Nominal pump flows [cm^3/s]
Fbar = [175;175];        % Mean disturbance flows [cm^3/s]

options = optimoptions('fsolve','Display','none');
xs = fsolve(@(x) FourTankSystem(0,x,us,Fbar,p),10000*ones(4,1),options);
hs = FourTankSensors(xs,p);

%% Problem 12.1 - Simulate the continuous-discrete plant

Ts = 4;                 % Sensor sampling period [s]
dt = 0.1;               % Euler-Maruyama integration step [s]
N = 300;
T = 0:Ts:N*Ts;

tauF = 40;
aF = 1/tauF;
stdF = 17.5;             % Stationary OU standard deviation [cm^3/s]
sigmaOU = stdF*sqrt(2*aF);
R = 3*eye(4);           % Measurement covariance [cm^2]

% Apply two separate 25% pump steps and return to the nominal flows.
U = repmat(us,1,N);
tInput = T(1:N);
U(1,tInput >= 160 & tInput < 520) = 1.25*us(1);
U(2,tInput >= 600 & tInput < 1000) = 1.25*us(2);

rng(1,'twister')
[~,XA,T,Y,Ztrue] = SimulateFourTankSDE(xs,U,Fbar,R,p,Ts,dt,'ou',aF,sigmaOU);

% Save the simulation truth at the sensor sampling instants.
Xtrue = XA(:,1:round(Ts/dt):end);
Htrue = [Ztrue; Xtrue(3:4,:)./(rho*A_tank(3:4))];
Ftrue = Xtrue(5:6,:);

%% Problem 12.2 - Continuous-discrete extended Kalman filter

% Perturb the initial estimate to show how the filter recovers.
xhat = [rho*A_tank.*(hs + [4;-4;3;-3]); Fbar + [10;-10]];
P = diag([5*rho*A_tank; stdF*ones(2,1)].^2);

Xhat = zeros(6,N+1);
E = zeros(4,N+1);
stdState = zeros(6,N+1);
stdInnovation = zeros(4,N+1);

for k = 1:N+1
    if k == 1
        uPrevious = us;
        elapsedTime = 0;
    else
        uPrevious = U(:,k-1);
        elapsedTime = Ts;
    end

    [xhat,P,e,Re] = FourTankCDEKF(Y(:,k),uPrevious,xhat,P,p,Fbar,aF,sigmaOU,elapsedTime,R);

    Xhat(:,k) = xhat;
    E(:,k) = e;
    stdState(:,k) = sqrt(diag(P));
    stdInnovation(:,k) = sqrt(diag(Re));
end

Hhat = Xhat(1:4,:)./(rho*A_tank);
Fhat = Xhat(5:6,:);
errorState = [Hhat-Htrue; Fhat-Ftrue];
stdError = [stdState(1:4,:)./(rho*A_tank); stdState(5:6,:)];

%% Simulation results

rmseSensor = sqrt(mean((Y-Htrue).^2,2));
rmseEstimate = sqrt(mean(errorState.^2,2));
rmseNominal = sqrt(mean((Fbar-Ftrue).^2,2));
meanInnovation = mean(E,2);
insideInnovationBand = 100*mean(abs(E) <= 2*stdInnovation,2);

fprintf('OU diffusion coefficient [cm^3/s^(3/2)]: %.4f\n',sigmaOU)

Tank = (1:4)';
disp(table(Tank,rmseSensor,rmseEstimate(1:4),meanInnovation, ...
    insideInnovationBand,'VariableNames',{'Tank','SensorRMSE_cm', ...
    'EstimateRMSE_cm','MeanInnovation_cm','WithinTwoSigma_pct'}))

Disturbance = (3:4)';
disp(table(Disturbance,rmseNominal,rmseEstimate(5:6), ...
    'VariableNames',{'Flow','NominalRMSE_cm3_per_s','EstimateRMSE_cm3_per_s'}))

%% Figures

figureDir = fullfile(fileparts(mfilename('fullpath')),'figures');
if ~isfolder(figureDir)
    mkdir(figureDir)
end

stateName = {'h_1','h_2','h_3','h_4','F_3','F_4'};
stateUnit = {'cm','cm','cm','cm','cm^3/s','cm^3/s'};
truth = [Htrue;Ftrue];
estimate = [Hhat;Fhat];
minutes = T/60;
estimateColor = [0.00 0.30 0.65];
bandColor = [0.63 0.27 0.00];

fig1 = figure('Color','w','Name','True and estimated states', ...
    'Position',[80 80 1060 660]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')
for i = 1:6
    nexttile
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',12)
    hold on
    if i <= 4
        plot(minutes,Y(i,:),'.','Color',[0.72 0.72 0.72], ...
            'MarkerSize',5,'DisplayName','Measurement')
    end
    plot(minutes,truth(i,:),'k-','LineWidth',1.25,'DisplayName','True')
    plot(minutes,estimate(i,:),'--','Color',estimateColor, ...
        'LineWidth',1.4,'DisplayName','CD-EKF')
    ylabel(sprintf('%s [%s]',stateName{i},stateUnit{i}))
    grid on
    if i >= 5
        xlabel('Time [min]')
    end
    if i == 1
        legend('Location','best','Color','w','TextColor','k')
    end
end
exportgraphics(fig1,fullfile(figureDir,'problem12_estimates.png'),'Resolution',180)

fig2 = figure('Color','w','Name','Estimation errors', ...
    'Position',[80 80 1060 660]);
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')
for i = 1:6
    nexttile
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',12)
    hold on
    yline(0,':','Color',[0.35 0.35 0.35],'HandleVisibility','off')
    plot(minutes,errorState(i,:),'-','Color',estimateColor, ...
        'LineWidth',1.2,'DisplayName','Estimation error')
    plot(minutes,2*stdError(i,:),'--','Color',bandColor, ...
        'LineWidth',1.15,'DisplayName','Two standard deviations')
    plot(minutes,-2*stdError(i,:),'--','Color',bandColor, ...
        'LineWidth',1.15,'HandleVisibility','off')
    ylabel(sprintf('%s error [%s]',stateName{i},stateUnit{i}))
    grid on
    if i >= 5
        xlabel('Time [min]')
    end
    if i == 1
        legend('Location','best','Color','w','TextColor','k')
    end
end
exportgraphics(fig2,fullfile(figureDir,'problem12_errors.png'),'Resolution',180)

fig3 = figure('Color','w','Name','Measurement innovations', ...
    'Position',[80 80 1060 520]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
for i = 1:4
    nexttile
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',12)
    hold on
    yline(0,':','Color',[0.35 0.35 0.35],'HandleVisibility','off')
    plot(minutes,E(i,:),'-','Color',estimateColor, ...
        'LineWidth',1.0,'DisplayName','Innovation')
    plot(minutes,2*stdInnovation(i,:),'--','Color',bandColor, ...
        'LineWidth',1.15,'DisplayName','Two standard deviations')
    plot(minutes,-2*stdInnovation(i,:),'--','Color',bandColor, ...
        'LineWidth',1.15,'HandleVisibility','off')
    ylabel(sprintf('e_%d [cm]',i))
    grid on
    if i >= 3
        xlabel('Time [min]')
    end
    if i == 1
        legend('Location','best','Color','w','TextColor','k')
    end
end
exportgraphics(fig3,fullfile(figureDir,'problem12_innovations.png'),'Resolution',180)
