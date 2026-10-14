function [x_filt, P_filt, x_pred, P_pred, e, Re, Kfx] = kalman_filter_time_variant(y_k, u_k, x_pred_k, P_pred_k, A, B, C, D, G, Q, R, S)
% KALMAN_FILTER_TIME_VARIANT Discrete-Time Time-Variant Kalman Filter
%
% Syntax:
%   [x_filt, P_filt, x_pred, P_pred, e, Re, Kfx] = ...
%       kalman_filter_time_variant(y_k, u_k, x_pred_k, P_pred_k, A, B, C, D, G, Q, R, S)

if nargin < 12 || isempty(S)
    S = zeros(size(G, 2), size(R, 1));
end

% 1. Measurement Update
y_hat = C * x_pred_k + D * u_k;          
e     = y_k - y_hat;                      
Re    = C * P_pred_k * C' + R;            

if rcond(Re) < 1e-12
    error('Innovation covariance Re is singular or ill-conditioned.');
end

Kfx   = (P_pred_k * C') / Re;             
Kfw   = S / Re;                           

x_filt = x_pred_k + Kfx * e;              
P_filt = P_pred_k - Kfx * Re * Kfx';  
P_filt = (P_filt + P_filt') / 2;

w_filt = Kfw * e;                         
Q_filt = Q - Kfw * Re * Kfw';

% 2. Time Update (One-step prediction)
x_pred = A * x_filt + B * u_k + G * w_filt;
P_pred = A * P_filt * A' + G * Q_filt * G' - A * Kfx * S' * G' - G * S * Kfx' * A';
P_pred = (P_pred + P_pred') / 2;          

end