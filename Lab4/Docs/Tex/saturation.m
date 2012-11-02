function Im = saturation(Vcmd, Vsupply, thetadot, Kb, Ka, Rm)
% determines motor current as a function of the command voltage, supply
% voltages, angular speed, and motor back emf constant
% Vsupply are the supplied votages to the power supply
% thetadot is the current angular speed of the motor
% Kb is the back emf constant
% Ka is the nominal motor driver gain
% R is nominanal motor resistance

% for thetadot no load = 600 rad/s then Kb = 0.04

% for the given command current, how much voltage would we need
Vb = Kb*thetadot;
Icmd = Ka * Vcmd;
Vrequired = Vb + Rm * Icmd;

if abs(Vrequired) < Vsupply
    Im = Icmd; % if we can supply enough voltage, then we can send this much current
else
    % otherwise, need to figure out what current we are limited to
    if Vcmd > 0
        Im = (Vsupply - Vb) / Rm;
    else
        Im = (-Vsupply - Vb) / Rm;
    end
end
end