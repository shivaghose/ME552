function [F] = estimateForce(I, x)
% numerically computes the partial derivative of flux with respect to x
% to estimate the force exerted

% expect values in the range of 0.1 N

N = 1000; % number of turns on the electromagnet TODO UNKNOWN

dx = 1e-6;

[phimain0 phileak0] = estimateFlux(I, x);

[phimain1 phileak1] = estimateFlux(I, x + dx);

deltaphi = (phimain1 + phileak1) - (phimain0 + phileak0);

dphidx = deltaphi/dx;

F = 0.5 * N * I * dphidx;

end

