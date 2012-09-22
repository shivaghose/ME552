function [avgC] = estimateMagnetParametersDumb()
% estimates C for the really dumb model

% actual data
%X = [0.52; 1.95; 2.38; 2.88; 3.38; 3.88; 4.31] ./ 1000; % want m
%V = [0.74; 1.17; 1.349; 1.565; 1.74; 2; 2.25];

% compute current
%I = driverVtoI(V);

% different data set
I = [0.358, 0.308, 0.365, 0.298, 0.288, 0.277] ;
X = [3.23, 2.73, 3.235, 3.13, 3.035, 2.935] ./ 1000;

C = zeros(length(X),1);

f = -0.01 * 9.81; % 10 grams

for i = 1:length(C)
    C(i) = -(f * X(i) * X(i)) / (I(i) * I(i));
end
C

avgC = mean(C);
end

