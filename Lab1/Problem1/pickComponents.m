function [Rs] = pickComponents(Clead, Clag, zlead, plead, zlag, plag, Gain, R1amp)
% where Clead and Clag are the known capacitances, z and p are the target
% zero and pole locations, and gain is the target gain 
% note that for this problem want
% zlead = 10; plead = 100; zlag = 0.2; plag = 0.1; Gain = 5;
% note that you give it the negative of the poles


% compute resistances for lead filter
R1lead = 1 / (Clead * zlead);
R2lead = 1 / ((plead - zlead) * Clead);

% compute resistances for lag filter
R2lag = 1 / (Clag * zlag);
R1lag = 1 / (Clag * plag) - R2lag;

% compute the gain of the system
K = R2lag / (R1lag + R2lag);

Ratio = Gain / K; % this is how much amplification we need to get the target gain

R2amp = (Ratio - 1)*R1amp;

Rs = [R1lead; R2lead; R1lag; R2lag; R2amp];
end

