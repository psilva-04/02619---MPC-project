function [T,X,Y,Z,D] = SimulateFourTankPWC(x0,U,ds,Qd,R,p,Ts)

N = size(U,2);

T = 0:Ts:N*Ts;

X = zeros(4,N+1);
Y = zeros(4,N+1);
Z = zeros(2,N+1);
D = zeros(2,N);

X(:,1) = x0;

Ld = chol(Qd,'lower');
Lv = chol(R,'lower');

Y(:,1) = FourTankSensors(X(:,1),p) + Lv*randn(4,1);
Z(:,1) = FourTankOutputs(X(:,1),p);

holdSteps = 10;

for k = 1:N

    u = U(:,k);

    if k == 1 || mod(k-1,holdSteps) == 0
        d = ds + Ld*randn(2,1);
    end

    D(:,k) = d;

    [~,x] = ode45( ...
        @(t,x) FourTankSystem(t,x,u,d,p), ...
        [T(k) T(k+1)], ...
        X(:,k));

    X(:,k+1) = x(end,:)';

    Y(:,k+1) = FourTankSensors(X(:,k+1),p) ...
             + Lv*randn(4,1);

    Z(:,k+1) = FourTankOutputs(X(:,k+1),p);

end

end