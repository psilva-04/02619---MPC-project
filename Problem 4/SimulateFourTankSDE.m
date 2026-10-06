function [t,XA,T,Y,Z] = SimulateFourTankSDE( ...
    x0,U,Fbar,RY, RZ,p,Ts,dt,model,aF,sigmaF)
% Simulates the stochastic nonlinear four-tank model.
%
% The augmented state is (as shown in fourtankSDE)
%
% XA = [m1; m2; m3; m4; F3; F4]
%
% The continuous-time model
%
% dXA = f dt + sigma dW
%
% is simulated using Euler-Maruyama.

%% Simulation dimensions

N = size(U,2);

stepsPerSample = round(Ts/dt);
nSteps = N*stepsPerSample;

t = 0:dt:N*Ts;
T = 0:Ts:N*Ts;


%% Augmented state

XA = zeros(6,nSteps+1);

XA(:,1) = [
    x0
    Fbar
];


%% Measurements and controlled outputs

Y = zeros(4,N+1);
Z = zeros(2,N+1);

LvY = chol(RY,'lower');
LvZ = chol(RZ,'lower');


%% Initial measurement

x = XA(1:4,1);

Y(:,1) = FourTankSensors(x,p) ...
       + LvY*randn(4,1);

Z(:,1) = FourTankOutputs(x,p) ...
        + LvZ*randn(2,1);


%% Euler-Maruyama simulation

for k = 1:N

    % Piecewise-constant manipulated input
    u = U(:,k);


    for j = 1:stepsPerSample

        i = (k-1)*stepsPerSample + j;

        xa = XA(:,i);


        %% Continuous-time SDE

        [f,sigma] = FourTankSDE( ...
            t(i), ...
            xa, ...
            u, ...
            Fbar, ...
            aF, ...
            sigmaF, ...
            p, ...
            model);


        %% Wiener increment

        % dW ~ N(0, I*dt)

        dW = sqrt(dt)*randn(2,1);


        %% Euler-Maruyama

        % xa(k+1) = xa(k) + f*dt + sigma*dW

        XA(:,i+1) = xa ...
                  + f*dt ...
                  + sigma*dW;

    end


    %% Measurement at sampling instant

    iSample = k*stepsPerSample + 1;

    x = XA(1:4,iSample);

    Y(:,k+1) = FourTankSensors(x,p) ...
             + LvY*randn(4,1);

    Z(:,k+1) = FourTankOutputs(x,p) ...
             + LvZ*randn(2,1);

end

end