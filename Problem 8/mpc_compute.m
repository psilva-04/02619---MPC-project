function [uk, Uk, Zk, info] = mpc_compute(mpc, x_k, w_k, u_prev, R_k, U_bar, xinit)
    
    % Free response
    bk = mpc.Phi_x * x_k + mpc.Phi_w * w_k;
    ck = R_k - bk;
    
    % Dynamic gradient vector g
    gz  = mpc.Mz * ck;
    gu  = mpc.Mu * U_bar;
    gdu = mpc.Mdu * u_prev;
    g = gz + gu + gdu;
    
    H     = mpc.H;
    l     = mpc.Umin;
    u     = mpc.Umax;
    A_qp  = mpc.Lambda;
    bl_qp = mpc.dUmin + mpc.I0 * u_prev;
    bu_qp = mpc.dUmax + mpc.I0 * u_prev;

    [Uk, info] = qpsolver(H, g, l, u, A_qp, bl_qp, bu_qp, xinit);
    
    
    uk = Uk(1:mpc.nu);
    Zk = bk + mpc.Gamma * Uk;
    
end