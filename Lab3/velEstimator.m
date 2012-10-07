function [thetadot] = velEstimator(T,Theta,I,J,Kt,B,Ff,thetaZ,thetadotZ,estThetaDot)
% estimate current Velocity?

% need to integrate current
SumI = 0;
for i = 1:(length(T) - 1)
    SumI = SumI + I(i) * (T(i + 1) - T(i));
end

thetadot = (Kt / J) * SumI - (B / J) * Theta(end) + (B / J) * thetaZ + J * thetadotZ + integrateCoulombicFriction(T,estThetaDot,I,Kt,Ff);

end

function [SumF] = integrateCoulombicFriction(T,estThetaDot,I,Kt,Ff)
SumF = 0;

for i = 1:(length(T) - 2) % don't have estThetaDot for the current time step
    friction = coulombicFriction(Ff,Kt*I(i),estThetaDot(i));
    SumF = SumF + friction * (T(i + 1) - T(i));
end

% repeat the last time step
if (length(T) > 1)
    friction = coulombicFriction(Ff,Kt*I(end),estThetaDot(end));
    SumF = SumF + friction * (T(end) - T(end-1));
end

end

function [torque] = coulombicFriction(Ff, motorTorque, thetadot)
% coulombic friction cannont exceed the magnitude of motor torque and
% opposes attempted motion

% simplified friction model which uses the same constant for kinetic and static friction

if (abs(thetadot) < 1e-6) % if not moving, friction opposes attempted motion
    if abs(motorTorque) < Ff
        FF = abs(motorTorque); % and doesn't exceed the magnitude of the applied torque
    else
        FF = Ff;
    end
    torque = - FF * sign(motorTorque);
    
else % if we are moving oppose the actual motion with kinetic friciton
    torque = - Ff * sign(thetadot);
end

end