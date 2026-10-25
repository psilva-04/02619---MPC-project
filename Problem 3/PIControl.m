function [u, i] = PIControl(i, r, y, us, Kc, Ti, Ts)
    e = r - y;
    i = i + (Kc*Ts/Ti) * e;
    u = us + Kc*e + i;
end