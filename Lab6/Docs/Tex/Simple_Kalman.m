function [ x_new, p_new ] = Simple_Kalman(measurement, x_old, p, r, q)
%SIMPLE_KALMAN Simple Kalman filter for a 1-D signal.
%   The Kalman filter uses an uncertainty estimate of the system noise to
%   regulate the signal and mitigate the effects of noise. The Kalman
%   update is as follows:
%       1) Update Kalman gain (K)
%       2) Update the 1-D State Variable (x)
%       3) Update the estimated error (p)
%   The system requires an estimate of the System Error and uses both the 
%   process noise (q) and the sensor noise (r) to determine the new state
%   (x_new). Run the Simple_Kalman_Init script to set the initial values of
%   the state vairbale (x) and the Estimated Error (p).

p_new = p + q;
K = p_new/(p_new+r);
x_new = x_old + K*(measurement - x_old);
p_new = (1-K)*p_new;
end