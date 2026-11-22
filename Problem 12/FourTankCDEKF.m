function [xhat,P,e,Re] = FourTankCDEKF(y,u,xhat,P,p,Fbar,aF,sigmaOU,Ts,R)
% Predict over the previous interval, then update using y at its end.
% At the first sample, use Ts = 0 to update the initial prior directly.

if Ts > 0
    [~,w] = ode45(@(t,w) predictionODE(t,w,u,p,Fbar,aF,sigmaOU), ...
        [0 Ts],[xhat;P(:)]);
    xhat = w(end,1:6)';
    P = reshape(w(end,7:end),6,6);
end
P = (P + P')/2;

% Measurement update (Lecture 6)
C = [diag(1./(p(12)*p(5:8))), zeros(4,2)];
yhat = FourTankSensors(xhat(1:4),p);
e = y - yhat;
Re = C*P*C' + R;
K = (P*C')/Re;

xhat = xhat + K*e;
P = P - K*Re*K';
P = (P + P')/2;
end

function wdot = predictionODE(t,w,u,p,Fbar,aF,sigmaOU)
xhat = w(1:6);
P = reshape(w(7:end),6,6);
[f,G] = FourTankSDE(t,xhat,u,Fbar,aF,sigmaOU,p,'ou');

% Use the Problem 5 Jacobian at the current estimated heights.
h = FourTankSensors(xhat(1:4),p);
lambda = p(1:4)./p(5:8).*sqrt(p(11)./(2*h));

A4 = [
    -lambda(1), 0,          lambda(3),  0
     0,        -lambda(2),  0,          lambda(4)
     0,         0,         -lambda(3),  0
     0,         0,          0,         -lambda(4)
];

Bd = p(12)*[zeros(2,2);eye(2)];
A = [A4,Bd;zeros(2,4),-aF*eye(2)];

Pdot = A*P + P*A' + G*G';
wdot = [f;Pdot(:)];
end
