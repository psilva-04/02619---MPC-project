function [uk, Uk, Zk, info] = mpc_compute(mpc, x_k, w_k, u_prev, R_k, U_bar, xinit)
    
    % Free response
    bk = mpc.Phi_x * x_k + mpc.Phi_w * w_k;
    ck = R_k - bk;
    
    % Dynamic gradient vector g
    gz  = mpc.Mz * ck;
    gu  = mpc.Mu * U_bar;
    gdu = mpc.Mdu * u_prev;
    g = gz + gu + gdu;

    gs = zeros(mpc.Nz, 1);
    gt = zeros(mpc.Nz, 1);
    g_bar = [g; gs; gt];
    
    H_qp  = mpc.H_bar;
    A_qp  = mpc.A_bar;
    l_qp  = mpc.l_bar;
    u_qp  = mpc.u_bar;

    Rmin_k = repmat(mpc.zmin, mpc.N, 1); 
    Rmax_k = repmat(mpc.zmax, mpc.N, 1);

    bl_qp = [mpc.dUmin + mpc.I0 * u_prev; Rmin_k - bk; -inf(mpc.Nz, 1)];
    bu_qp = [mpc.dUmax + mpc.I0 * u_prev; inf(mpc.Nz, 1); Rmax_k - bk];

    xinit_bar = [xinit; zeros(2 * mpc.Nz, 1)];

    [Uk_bar, info] = qpsolver(H_qp, g_bar, l_qp, u_qp, A_qp, bl_qp, bu_qp, xinit_bar);
    
    Uk = Uk_bar(1:mpc.N * mpc.nu);
    uk = Uk(1:mpc.nu);
    Zk = bk + mpc.Gamma * Uk;
end