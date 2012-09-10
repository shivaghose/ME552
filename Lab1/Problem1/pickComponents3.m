function [Rs] = pickComponents3(C1, C2)
% where Clead and Clag are the known capacitances, z and p are the target
% zero and pole locations, and gain is the target gain 
% note that for this problem want
% zlead = 10; plead = 100; zlag = 0.2; plag = 0.1; Gain = 5;
% note that you give it the negative of the poles

% revised yet again for another circuit
% hopefully this actually works!

% parameters for our desired system
z1 = 10; p1 = 100; z2 = 0.2; p2 = 0.1; Gain = 5;

% lead component
R2 = 1 / (C2 * z1);
R1 = 1/ (C1 * p1);

% lag component
R3 = 1 / (z2 * C1) - R1;
R4 = 1 / (p2 * C2) - R2;

curGain = - (R4 * R2 * C2 * (R1 * C1 + R3 * C1)) / (R3 * R1 * C1 * (R4 * C2 + R2 * C2))

ratio = Gain / curGain

Rs = [R1; R2; R3; R4];

% try flipping it the other way doesn't work, requires negative resistances
% % lag component
% R2 = 1 / (C2 * z2);
% R1 = 1/ (C1 * p2);
% 
% % lead component
% R3 = 1 / (z1 * C1) - R1;
% R4 = 1 / (p1 * C2) - R2;
% 
% curGain = - (R4 * R2 * C2 * (R1 * C1 + R3 * C1)) / (R3 * R1 * C1 * (R4 * C2 + R2 * C2))
% 
% ratio = Gain / curGain
% 
% Rs = [R1; R2; R3; R4]
end

