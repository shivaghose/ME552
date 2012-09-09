function [Vout, I, P] = simulate(R, V)
% Runs a simulation of the circuit to calculate currents and power
% dissipated through each resistor in the circuit

% R = [R1, R2, R3, R4, Rf, Rg] expects the full R
% V = [V1, V2, V3, V4]
% extract values
    R1 = R(1); R2 = R(2); R3 = R(3); R4 = R(4); Rf = R(5); Rg = R(6);
    V1 = V(1); V2 = V(2); V3 = V(3); V4 = V(4);

% compute Vout
a = computeA(R);
Vout = a(1) * V(1) + a(2) * V(2) - a(3) * V(3) - a(4) * V(4);

% compute Vplus
Vplus = computeVplus(R,V);

% now compute currents and power
I1 = (V1 - Vplus) / R1; % resistor 1
P1 = I1 * (V1 - Vplus);

I2 = (V2 - Vplus) / R2; % resistor 2
P2 = I2 * (V2 - Vplus);

Ig = Vplus / Rg; % resistor to ground
Pg = Ig * Vplus;

I3 = (V3 - Vplus) / R3; % resistor 3
P3 = I3 * (V3 - Vplus);

I4 = (V4 - Vplus) / R4; % resistor 4
P4 = I4 * (V4 - Vplus);

If = (Vplus - Vout) / Rf;
Pf = If * (Vplus - Vout);

I = [I1; I2; I3; I4; If; Ig];
P = [P1; P2; P3; P4; Pf; Pg];
end

% computes voltage at non-inverting terminal of op-amp
function [Vplus] = computeVplus(R,V)
    % extract values
    R1 = R(1); R2 = R(2); R3 = R(3); R4 = R(4); Rf = R(5); Rg = R(6);
    V1 = V(1); V2 = V(2);
    
    Vplus = (1 / (R1 / Rg + R1 / R2 + 1)) * V1 + (1 / (R2 / R1 + R2 / Rg + 1)) * V2;
end

% for computing a
function [a] = computeA(R) 
    a = [computeA1(R); computeA2(R); computeA3(R); computeA4(R)];
end

function [a1] = computeA1(R) 
    % extract values
    R1 = R(1); R2 = R(2); R3 = R(3); R4 = R(4); Rf = R(5); Rg = R(6);
    
    % compute a1
    a1 = (Rf / R3 + Rf / R4 + 1) * 1 / (R1 / Rg + R1 / R2 + 1);
end

function [a2] = computeA2(R) 
    % extract values
    R1 = R(1); R2 = R(2); R3 = R(3); R4 = R(4); Rf = R(5); Rg = R(6);
    
    a2 = (Rf / R3 + Rf / R4 + 1) * 1 / (R2 / R1 + R2 / Rg + 1);
    
end

function [a3] = computeA3(R)
    % extract values
    R3 = R(3); Rf = R(5);
    
    a3 = Rf / R3;
end

function [a4] = computeA4(R)
    % extract values
    R4 = R(4); Rf = R(5);
    
    a4 = Rf / R4;
end