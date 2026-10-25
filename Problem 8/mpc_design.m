function mpc = mpc_design(Ad, Bd, Gd, Cz, N, Wz, Wu, Wdu, umin, umax, dumin, dumax)

    nx = size(Ad, 1);
    nu = size(Bd, 2);
    nz = size(Cz, 1);
    
    % Build prediction matrices
    [Phi_x, Phi_w, Gamma] = build_prediction_matrices(Ad, Bd, Gd, Cz, N);

    % Weights for Hessian
    I_N = eye(N);
    Wz_bar = kron(I_N, Wz);
    Wu_bar = kron(I_N, Wu);
    Wdu_bar = kron(I_N, Wdu);
    
    I_nu = eye(nu);
    Lambda = kron(eye(N), I_nu) - kron(diag(ones(N-1,1), -1), I_nu);
    I0 = zeros(N*nu, nu);
    I0(1:nu, 1:nu) = I_nu;

    % Hessian calculation (only do once)
    Hz = (Wz_bar*Gamma)'*(Wz_bar*Gamma);
    Hu = Wu_bar'*Wu_bar;
    Hdu = (Wdu_bar*Lambda)'*(Wdu_bar*Lambda);
    H = Hz + Hu + Hdu;
    
    % Gradient multiplication matrices
    Mz  = -(Wz_bar * Gamma)' * Wz_bar;
    Mu  = -Wu_bar' * Wu_bar;
    Mdu = -(Wdu_bar * Lambda)' * Wdu_bar * I0;
    
    mpc.Phi_x  = Phi_x;
    mpc.Phi_w  = Phi_w;
    mpc.Gamma  = Gamma;
    mpc.H      = H;
    mpc.Mz     = Mz;
    mpc.Mu     = Mu;
    mpc.Mdu    = Mdu;
    mpc.Lambda = Lambda;
    mpc.I0     = I0;
    mpc.Umin   = repmat(umin(:),  N, 1);
    mpc.Umax   = repmat(umax(:),  N, 1);
    mpc.dUmin  = repmat(dumin(:), N, 1);
    mpc.dUmax  = repmat(dumax(:), N, 1);
    mpc.nu     = nu;
    mpc.nz     = nz;
    mpc.N      = N;

end

function [Phi_x, Phi_w, Gamma] = build_prediction_matrices(Ad, Bd, Gd, Cz, N)

    nx = size(Ad, 1);
    nu = size(Bd, 2);
    nw = size(Gd,2);
    nz = size(Cz, 1);
    
    Phi_x = zeros(N * nz, nx);
    Phi_w = zeros(N * nz, nw);
    Gamma = zeros(N * nz, N * nu);
    
    A_pow = eye(nx);
    
    for i = 1:N
    
        row_idx = (i-1)*nz + 1 : i*nz;
    
        Hi    = Cz * A_pow * Bd;
        Phi_w(row_idx, :) = Cz * A_pow * Gd; % C*G, C*A*G, C*A^2*G...
    
        A_pow = A_pow * Ad;
        Phi_x(row_idx, :) = Cz * A_pow;
    
        for j = 1:(N - i + 1)
            row_block = (i + j - 2)*nz + 1 : (i + j - 1)*nz;
            col_block = (j - 1)*nu + 1 : j*nu;
            Gamma(row_block, col_block) = Hi;
        end
    end
end
