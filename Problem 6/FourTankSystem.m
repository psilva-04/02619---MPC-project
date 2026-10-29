function dx = FourTankSystem(t,x,u,d,p)
% Deterministic nonlinear dynamics of the modified four-tank system.
% p = [a1 a2 a3 a4 A1 A2 A3 A4 gamma1 gamma2 g rho]'

a     = p(1:4);
A     = p(5:8);
gamma = p(9:10);
g     = p(11);
rho   = p(12);

F1 = u(1);
F2 = u(2);

F3 = d(1);
F4 = d(2);

h = x./(rho*A);

qout = a.*sqrt(2*g*h);

q1in = gamma(1)*F1;
q2in = gamma(2)*F2;
q3in = (1-gamma(2))*F2;
q4in = (1-gamma(1))*F1;

dx = rho * [
    q1in + qout(3) - qout(1)
    q2in + qout(4) - qout(2)
    q3in + F3      - qout(3)
    q4in + F4      - qout(4)
];

end