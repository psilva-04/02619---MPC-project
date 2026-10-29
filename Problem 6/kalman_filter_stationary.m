function [x_filt, x_pred, e] = kalman_filter_stationary(y_k, u_k, x_pred_k, A, B, C, D, G, Kfx, Kfw)
% KALMAN_FILTER_STATIONARY Discrete-Time Stationary Kalman Filter
%
% Syntax:
%   [x_filt, x_pred, e] = ...
%       kalman_filter_stationary(y_k, u_k, x_pred_k, A, B, C, D, G, Kfx, Kfw)

if nargin < 10 || isempty(Kfw)
    Kfw = zeros(size(G, 2), length(y_k));
end

% 1. Innovation
y_hat = C * x_pred_k + D * u_k;
e     = y_k - y_hat;

% 2. Filtered State and Process Noise
x_filt = x_pred_k + Kfx * e;
w_filt = Kfw * e;

% 3. State Prediction to k+1|k
x_pred = A * x_filt + B * u_k + G * w_filt;

end