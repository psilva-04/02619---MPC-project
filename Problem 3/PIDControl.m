function [u, i, e] = PIDControl(i, e_prev, r, y, us, Kc, Ti, Td, Ts)
    
    e = r - y;
    P = Kc * e; 
    i = i + (Kc * Ts / Ti) * e; 
    D = Kc * Td * (e - e_prev) / Ts;
    u = us + P + i + D;

end