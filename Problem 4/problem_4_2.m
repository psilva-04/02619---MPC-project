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

zs = FourTankOutputs(xs,p);

% Simulation settings

tspan = [0 1000];

steps = [0.10 0.25 0.50];

nSteps = length(steps);

%% 4.2 - Stochastic nonlinear step responses with OU disturbances

Ts = 4;
dt = 0.1;

N = 250;        % 250*4 = 1000 s
T = 0:Ts:N*Ts;

t_step = 52;    % Must coincide with Ts = 4 sampling grid

tauF = 40;
aF = 1/tauF;


% Noise levels

processNoise = [0.02 0.05 0.10];

sensorNoise = [1 2 3];

noiseNames = {
    'Low'
    'Medium'
    'High'
};


% STOCHASTIC SIMULATIONS
% Low noise:    Process = 2%,  Sensor R = I
% Medium noise: Process = 5%,  Sensor R = 2I
% High noise:   Process = 10%, Sensor R = 3I

Gstoch1 = cell(2,2,nSteps,3);
Gstoch2 = cell(2,2,nSteps,3);
DataStoch = cell(2,2,nSteps,3);

for n = 1:3

    % Process noise

    stdF = processNoise(n)*175;
    sigmaOU = stdF*sqrt(2*aF);


    % Sensor noise

    RY = sensorNoise(n)*eye(4);
    RZ = sensorNoise(n)*eye(2);


    % Step Responses in F1

    figure('Name',sprintf( ...
    '%s Noise - F1 Steps', ...
    noiseNames{n}), ...
    'NumberTitle','off')

    ZnormF1 = zeros(2,N+1,nSteps);

    for i = 1:nSteps
    
        stepSize = steps(i);
    
        % Nominal input for entire simulation
        U = repmat(us,1,N);
    
        % Apply step in F1 at t = 52 s
        U(1,T(1:N) >= t_step) = us(1)*(1 + stepSize);
    
        rng(1)
    
        [tOU,XAOU,TOU,YOU,ZOU] = ...
            SimulateFourTankSDE( ...
            xs, ...
            U, ...
            Fbar, ...
            RY, ...
            RZ, ...
            p, ...
            Ts, ...
            dt, ...
            'ou', ...
            aF, ...
            sigmaOU);
    
        delta_u = us(1)*stepSize;
        ZnormF1(:,:,i) = (ZOU - zs)/delta_u;

        subplot(3,1,i)
    
        plot(TOU,ZOU,'LineWidth',1.3)
        grid on
    
        xlabel('Time [s]')
        ylabel('Output [cm]')
    
        title(sprintf('F1 step = +%.0f%%',100*stepSize))
    
        legend('h_1','h_2', ...
            'Location','best')
    
        upper_limit = (floor(max(ZOU(:))/10) + 1)*10;
        ylim([60 upper_limit])
    
        xline(t_step,'--','Step','HandleVisibility','off')
    
    end
    
    
    sgtitle(sprintf( ...
    'F1 Steps - %s Noise: Process = %.0f%%, Sensor R = %d%%', ...
    noiseNames{n}, ...
    100*processNoise(n), ...
    sensorNoise(n)))

    figure('Name',sprintf( ...
    '%s Noise - Normalized F1 Steps', ...
    noiseNames{n}), ...
    'NumberTitle','off')

    subplot(2,1,1)
    hold on
    
    for i = 1:nSteps
        plot(TOU,ZnormF1(1,:,i),'LineWidth',1.3)
    end
    
    grid on
    xlabel('Time [s]')
    ylabel('\Deltah_1 / \DeltaF_1')
    legend('10%','25%','50%','Location','best')
    title('Normalized response of h_1')
    
    subplot(2,1,2)
    hold on
    
    for i = 1:nSteps
        plot(TOU,ZnormF1(2,:,i),'LineWidth',1.3)
    end
    
    grid on
    xlabel('Time [s]')
    ylabel('\Deltah_2 / \DeltaF_1')
    legend('10%','25%','50%','Location','best')
    title('Normalized response of h_2')
    
    sgtitle(sprintf( ...
        'Normalized F1 Steps - %s Noise', ...
        noiseNames{n}))


    % Step Responses in F2
    
    figure('Name',sprintf( ...
    '%s Noise - F2 Steps', ...
    noiseNames{n}), ...
    'NumberTitle','off')
    
    ZnormF2 = zeros(2,N+1,nSteps);

    for i = 1:nSteps
    
        stepSize = steps(i);
    
        % Nominal input for entire simulation
        U = repmat(us,1,N);
    
        % Apply step in F2 at t = 52 s
        U(2,T(1:N) >= t_step) = us(2)*(1 + stepSize);
    
        rng(1)
    
        [tOU,XAOU,TOU,YOU,ZOU] = ...
            SimulateFourTankSDE( ...
            xs, ...
            U, ...
            Fbar, ...
            RY, ...
            RZ, ...
            p, ...
            Ts, ...
            dt, ...
            'ou', ...
            aF, ...
            sigmaOU);
    
        delta_u = us(2)*stepSize;
        ZnormF2(:,:,i) = (ZOU - zs)/delta_u;

        subplot(3,1,i)
    
        plot(TOU,ZOU,'LineWidth',1.3)
        grid on
    
        xlabel('Time [s]')
        ylabel('Output [cm]')
    
        title(sprintf('F2 step = +%.0f%%',100*stepSize))
    
        legend('h_1','h_2', ...
            'Location','best')
    
        upper_limit = (floor(max(ZOU(:))/10) + 1)*10;
        ylim([60 upper_limit])
    
        xline(t_step,'--','Step','HandleVisibility','off')
    
    end
    
    
    sgtitle(sprintf( ...
    'F2 Steps - %s Noise: Process = %.0f%%, Sensor R = %d%%', ...
    noiseNames{n}, ...
    100*processNoise(n), ...
    sensorNoise(n)))

    figure('Name',sprintf( ...
    '%s Noise - Normalized F2 Steps', ...
    noiseNames{n}), ...
    'NumberTitle','off')

    subplot(2,1,1)
    hold on
    
    for i = 1:nSteps
        plot(TOU,ZnormF2(1,:,i),'LineWidth',1.3)
    end
    
    grid on
    xlabel('Time [s]')
    ylabel('\Deltah_1 / \DeltaF_2')
    legend('10%','25%','50%','Location','best')
    title('Normalized response of h_1')
    
    subplot(2,1,2)
    hold on
    
    for i = 1:nSteps
        plot(TOU,ZnormF2(2,:,i),'LineWidth',1.3)
    end
    
    grid on
    xlabel('Time [s]')
    ylabel('\Deltah_2 / \DeltaF_2')
    legend('10%','25%','50%','Location','best')
    title('Normalized response of h_2')
    
    sgtitle(sprintf( ...
        'Normalized F2 Steps - %s Noise', ...
        noiseNames{n}))

    % Plot disturbances F3 and F4

    figure('Name',sprintf( ...
        '%s Noise - OU Disturbances', ...
        noiseNames{n}), ...
        'NumberTitle','off')
    
    plot(tOU,XAOU(5,:),'LineWidth',1.3)
    hold on
    
    plot(tOU,XAOU(6,:),'LineWidth',1.3)
    
    yline(Fbar(1),'--','Nominal', ...
        'HandleVisibility','off')
    
    grid on
    
    xlabel('Time [s]')
    ylabel('Flow rate [cm^3/s]')
    
    legend('F_3','F_4', ...
        'Location','best')
    
    title(sprintf( ...
        'OU Disturbances - %s Noise (%.0f%%)', ...
        noiseNames{n}, ...
        100*processNoise(n)))

    % Transfer function identification

    t_start_id = t_step - 20;

    idx = TOU >= t_start_id;

    t_id = TOU(idx)';
    u_id = double(t_id >= t_step);

    for i = 1:nSteps

        % F1 -> h1, h2

        y1 = ZnormF1(1,idx,i)';
        y2 = ZnormF1(2,idx,i)';

        data11 = iddata(y1,u_id,Ts);
        data21 = iddata(y2,u_id,Ts);

        DataStoch{1,1,i,n} = data11;
        DataStoch{2,1,i,n} = data21;

        Gstoch1{1,1,i,n} = tfest(data11,1,0);
        Gstoch1{2,1,i,n} = tfest(data21,1,0);

        Gstoch2{1,1,i,n} = tfest(data11,2,0);
        Gstoch2{2,1,i,n} = tfest(data21,2,0);


        % F2 -> h1, h2

        y1 = ZnormF2(1,idx,i)';
        y2 = ZnormF2(2,idx,i)';

        data12 = iddata(y1,u_id,Ts);
        data22 = iddata(y2,u_id,Ts);

        DataStoch{1,2,i,n} = data12;
        DataStoch{2,2,i,n} = data22;

        Gstoch1{1,2,i,n} = tfest(data12,1,0);
        Gstoch1{2,2,i,n} = tfest(data22,1,0);

        Gstoch2{1,2,i,n} = tfest(data12,2,0);
        Gstoch2{2,2,i,n} = tfest(data22,2,0);

    end


    fprintf('\n%s noise\n',noiseNames{n})

    for i = 1:nSteps

        fprintf('\n%.0f%% step\n',100*steps(i))

        fprintf('\nG11: F1 -> h1\n')
        fprintf('First order:\n')
        disp(Gstoch1{1,1,i,n})
        fprintf('Second order:\n')
        disp(Gstoch2{1,1,i,n})

        fprintf('\nG12: F2 -> h1\n')
        fprintf('First order:\n')
        disp(Gstoch1{1,2,i,n})
        fprintf('Second order:\n')
        disp(Gstoch2{1,2,i,n})

        fprintf('\nG21: F1 -> h2\n')
        fprintf('First order:\n')
        disp(Gstoch1{2,1,i,n})
        fprintf('Second order:\n')
        disp(Gstoch2{2,1,i,n})

        fprintf('\nG22: F2 -> h2\n')
        fprintf('First order:\n')
        disp(Gstoch1{2,2,i,n})
        fprintf('Second order:\n')
        disp(Gstoch2{2,2,i,n})

    end


    for i = 1:nSteps

        figure('Name',sprintf( ...
            '%s Noise - TF Comparison - %.0f%% Step', ...
            noiseNames{n},100*steps(i)), ...
            'NumberTitle','off')

        subplot(2,2,1)
        compare( ...
            DataStoch{1,1,i,n}, ...
            Gstoch1{1,1,i,n}, ...
            Gstoch2{1,1,i,n})
        title('G_{11}: F_1 \rightarrow h_1')

        subplot(2,2,2)
        compare( ...
            DataStoch{1,2,i,n}, ...
            Gstoch1{1,2,i,n}, ...
            Gstoch2{1,2,i,n})
        title('G_{12}: F_2 \rightarrow h_1')

        subplot(2,2,3)
        compare( ...
            DataStoch{2,1,i,n}, ...
            Gstoch1{2,1,i,n}, ...
            Gstoch2{2,1,i,n})
        title('G_{21}: F_1 \rightarrow h_2')

        subplot(2,2,4)
        compare( ...
            DataStoch{2,2,i,n}, ...
            Gstoch1{2,2,i,n}, ...
            Gstoch2{2,2,i,n})
        title('G_{22}: F_2 \rightarrow h_2')

        sgtitle(sprintf( ...
            '%s Noise - First vs Second Order - %.0f%% Step', ...
            noiseNames{n},100*steps(i)))

    end
end