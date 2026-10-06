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
t_step = 50;

steps = [0.10 0.25 0.50];

nSteps = length(steps);

%% 4.1 - Deterministic nonlinear step responses (with constant disturbances)

% Step responses in F1

figure

ZnormF1 = cell(nSteps,1);
timeF1 = cell(nSteps,1);

for i = 1:nSteps

    stepSize = steps(i);
    
    u_after = us;
    u_after(1) = us(1)*(1 + stepSize);
    
    % Nominal input until t_step, then step
    ufun = @(t) us*(t < t_step) + u_after*(t >= t_step);
    
    [t,x] = ode45( ...
        @(t,x) FourTankSystem(t,x,ufun(t),Fbar,p), ...
        tspan, ...
        xs);

    % Calculate outputs
    z = zeros(length(t),length(zs));

    for k = 1:length(t)
        z(k,:) = FourTankOutputs(x(k,:)',p)';
    end

    delta_u = us(1)*stepSize;
    ZnormF1{i} = (z - zs')/delta_u;
    timeF1{i} = t;

    % Plot
    subplot(3,1,i)

    plot(t,z,'LineWidth',1.3)
    grid on

    xlabel('Time [s]')
    ylabel('Output [cm]')

    title(sprintf('F1 step = +%.0f%%',100*stepSize))

    legend('h_1','h_2', ...
        'Location','best')

    upper_limit = (floor(max(z(:))/10) + 1)*10;
    ylim([60 upper_limit])

    xline(t_step,'--','Step','HandleVisibility','off')

end

sgtitle('Deterministic Nonlinear Step Responses - Steps in F1')


% Normalized step responses in F1

figure

subplot(2,1,1)
hold on

for i = 1:nSteps
    plot(timeF1{i},ZnormF1{i}(:,1),'LineWidth',1.3)
end

grid on
xlabel('Time [s]')
ylabel('\Deltah_1 / \DeltaF_1')
legend('10%','25%','50%','Location','best')
title('Normalized response of h_1')

subplot(2,1,2)
hold on

for i = 1:nSteps
    plot(timeF1{i},ZnormF1{i}(:,2),'LineWidth',1.3)
end

grid on
xlabel('Time [s]')
ylabel('\Deltah_2 / \DeltaF_1')
legend('10%','25%','50%','Location','best')
title('Normalized response of h_2')

sgtitle('Normalized Deterministic Step Responses - Steps in F1')


% Step responses in F2

figure

ZnormF2 = cell(nSteps,1);
timeF2 = cell(nSteps,1);

for i = 1:nSteps

    stepSize = steps(i);

    u_after = us;
    u_after(2) = us(2)*(1 + stepSize);
    
    % Nominal input until t_step, then step
    ufun = @(t) us*(t < t_step) + u_after*(t >= t_step);
    
    [t,x] = ode45( ...
        @(t,x) FourTankSystem(t,x,ufun(t),Fbar,p), ...
        tspan, ...
        xs);

    % Calculate outputs
    z = zeros(length(t),length(zs));

    for k = 1:length(t)
        z(k,:) = FourTankOutputs(x(k,:)',p)';
    end

    delta_u = us(2)*stepSize;
    ZnormF2{i} = (z - zs')/delta_u;
    timeF2{i} = t;

    % Plot
    subplot(3,1,i)

    plot(t,z,'LineWidth',1.3)
    grid on

    xlabel('Time [s]')
    ylabel('Output [cm]')

    title(sprintf('F2 step = +%.0f%%',100*stepSize))

    legend('h_1','h_2', ...
        'Location','best')

    upper_limit = (floor(max(z(:))/10) + 1)*10;
    ylim([60 upper_limit])

    xline(t_step,'--','Step','HandleVisibility','off')

end

sgtitle('Deterministic Nonlinear Step Responses - Steps in F2')


% Normalized step responses in F2

figure

subplot(2,1,1)
hold on

for i = 1:nSteps
    plot(timeF2{i},ZnormF2{i}(:,1),'LineWidth',1.3)
end

grid on
xlabel('Time [s]')
ylabel('\Deltah_1 / \DeltaF_2')
legend('10%','25%','50%','Location','best')
title('Normalized response of h_1')

subplot(2,1,2)
hold on

for i = 1:nSteps
    plot(timeF2{i},ZnormF2{i}(:,2),'LineWidth',1.3)
end

grid on
xlabel('Time [s]')
ylabel('\Deltah_2 / \DeltaF_2')
legend('10%','25%','50%','Location','best')
title('Normalized response of h_2')

sgtitle('Normalized Deterministic Step Responses - Steps in F2')

%% Transfer function identification

Ts_id = 4;
t_start_id = t_step - 20;

Gdet1 = cell(2,2,nSteps);
Gdet2 = cell(2,2,nSteps);
DataDet = cell(2,2,nSteps);

for i = 1:nSteps

    t_id = (t_start_id:Ts_id:tspan(2))';
    u_id = double(t_id >= t_step);


    % F1 -> h1, h2

    t = timeF1{i};
    zNorm = ZnormF1{i};

    y1 = interp1(t,zNorm(:,1),t_id);
    y2 = interp1(t,zNorm(:,2),t_id);

    data11 = iddata(y1,u_id,Ts_id);
    data21 = iddata(y2,u_id,Ts_id);

    DataDet{1,1,i} = data11;
    DataDet{2,1,i} = data21;

    Gdet1{1,1,i} = tfest(data11,1,0);
    Gdet1{2,1,i} = tfest(data21,1,0);

    Gdet2{1,1,i} = tfest(data11,2,0);
    Gdet2{2,1,i} = tfest(data21,2,0);


    % F2 -> h1, h2

    t = timeF2{i};
    zNorm = ZnormF2{i};

    y1 = interp1(t,zNorm(:,1),t_id);
    y2 = interp1(t,zNorm(:,2),t_id);

    data12 = iddata(y1,u_id,Ts_id);
    data22 = iddata(y2,u_id,Ts_id);

    DataDet{1,2,i} = data12;
    DataDet{2,2,i} = data22;

    Gdet1{1,2,i} = tfest(data12,1,0);
    Gdet1{2,2,i} = tfest(data22,1,0);

    Gdet2{1,2,i} = tfest(data12,2,0);
    Gdet2{2,2,i} = tfest(data22,2,0);

end


for i = 1:nSteps

    fprintf('\n%.0f%% step\n',100*steps(i))

    fprintf('\nG11: F1 -> h1\n')
    fprintf('First order:\n')
    disp(Gdet1{1,1,i})
    fprintf('Second order:\n')
    disp(Gdet2{1,1,i})

    fprintf('\nG12: F2 -> h1\n')
    fprintf('First order:\n')
    disp(Gdet1{1,2,i})
    fprintf('Second order:\n')
    disp(Gdet2{1,2,i})

    fprintf('\nG21: F1 -> h2\n')
    fprintf('First order:\n')
    disp(Gdet1{2,1,i})
    fprintf('Second order:\n')
    disp(Gdet2{2,1,i})

    fprintf('\nG22: F2 -> h2\n')
    fprintf('First order:\n')
    disp(Gdet1{2,2,i})
    fprintf('Second order:\n')
    disp(Gdet2{2,2,i})

end


for i = 1:nSteps

    figure('Name',sprintf( ...
        'TF Comparison - %.0f%% Step', ...
        100*steps(i)), ...
        'NumberTitle','off')

    subplot(2,2,1)
    compare( ...
        DataDet{1,1,i}, ...
        Gdet1{1,1,i}, ...
        Gdet2{1,1,i})
    title('G_{11}: F_1 \rightarrow h_1')

    subplot(2,2,2)
    compare( ...
        DataDet{1,2,i}, ...
        Gdet1{1,2,i}, ...
        Gdet2{1,2,i})
    title('G_{12}: F_2 \rightarrow h_1')

    subplot(2,2,3)
    compare( ...
        DataDet{2,1,i}, ...
        Gdet1{2,1,i}, ...
        Gdet2{2,1,i})
    title('G_{21}: F_1 \rightarrow h_2')

    subplot(2,2,4)
    compare( ...
        DataDet{2,2,i}, ...
        Gdet1{2,2,i}, ...
        Gdet2{2,2,i})
    title('G_{22}: F_2 \rightarrow h_2')

    sgtitle(sprintf( ...
        'First vs Second Order Models - %.0f%% Step', ...
        100*steps(i)))

end