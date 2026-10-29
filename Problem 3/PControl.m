function u = PControl(r, y, us, Kc)
    e = r - y;
    u = us + Kc * e;
end

