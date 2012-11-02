function [alphadot, alphadotdot] = pendulumEQMotion(alpha, alphadot)
% Equation of motion for the pendulum by itself
% expects the following variables to be in the workspace
% g - gravity
% Iz2 - moment of inertia of pendulum about center of mass
% m2 - pendulum mass
% l2c - distance from axis of rotation to center of mass
% ba - viscous damping on pendulum
% Tca - torque exerted by static friction

sa = sin(alpha);

Tg = m2*g*l2c*sa;
Tc = coulombFriction(Tca, Tg, alphadot);

D = Iz2 + m2*l2c*l2c;

alphadotdot = (1/D)*(-ba*alphadot + m2*g*l2c*sa + Tc);
end

function [torque] = coulombFriction(Tf, externalTorque, thetadot)
% this friction model is making an assumption that isn't strictly true, we
% have static friction, but only the external torque exerted on each link
% will be considered.  This is not strictly correct, really we should
% consider the constraint forces as well
% it might not hold with the strength that it should
staticThreshold = 1e-6;
staticRatio = 1.1;

if (abs(thetadot) < staticThreshold) % if not moving, friction opposes attempted motion
    % if external torque does not exceed static friction
    if ((externalTorque < Tf) && (externalTorque >= 0)) || ((externalTorque > -Tf) && (externalTorque < 0))
        FF = abs(externalTorque);
    else
        FF = Tf;
    end
    torque = - staticRatio * FF * sign(externalTorque);
    
else % if we are moving oppose the actual motion with kinetic friciton
    if thetadot > staticThreshold
        torque = -Tf;
    else
        torque = Tf;
    end
end
end