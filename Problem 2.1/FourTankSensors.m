function y = FourTankSensors(x,p)
% Computes the measured liquid levels h1 and h2.
% p = [a1 a2 a3 a4 A1 A2 A3 A4 gamma1 gamma2 g rho]'

A   = p(5:8);
rho = p(12);

h = x./(rho*A);

y = [
    h(1)
    h(2)
];

end