function [] = manualParameterEntry(Tin,Thetain,Iin,Jin,Ktin,Bf,Br,Tff,Tfr,StaticFrictionMultiplier)
% we now have B and Ff from another data source in terms of Kt so they have
% been removed from our state description
% Let P be the vector of model parameter values
% P = [J]
% J: rotor moment of inertia
% B: viscous damping
% Fs: force of friction, not modeling static friction seperately
% Kt: motor torque constant

global Tobs Thetaobs Iobs J Kt Bforwards Tfforwards Breverse Tfreverse StaticFrictionX

Kt = Ktin; % now no longer a variable, single number!
Bforwards = Bf;
Breverse = Br;
Tfforwards = Tff;
Tfreverse = Tfr;
J = Jin;
StaticFrictionX = StaticFrictionMultiplier;

Tobs = Tin; Thetaobs = Thetain; Iobs = Iin; 

xstart = [0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial plot to make sure the shooting works correctly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T, X] = shootMotor(xstart);
Load = computeLoad(T,X);
figure
hold on;
plot(T,X(:,1),'b');
plot(Tobs,Thetaobs,'k');
plot(Tobs,Iobs,'r');
plot(T,Load,'g');
title('Initial Estimate curve');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
hold off;

end

% where x0 are the initial conditions initial angle, initial angular
% velocity, and P_ is the value of motor parameters to use
function [T, X] = shootMotor(xstart)

global Tobs

% determine starting and ending times from our data
tstart = 0;
tend = Tobs(end); % go quarter of a second longer past the end

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%[T, X] = ode45(@motorEQMotion, [tstart, tend], xstart, options);
T = linspace(tstart, tend, 10^4);
[X] = ode5(@motorEQMotion, T,xstart);

end

function [Load] = computeLoad(T,X)
global Bforwards Breverse
Load = zeros(length(T),1);
for i = 1:length(T)
    tau = motorTorque(motorCurrent(T(i)));
    thetadot = X(i,2);
    if (thetadot >= 0)
        B = Bforwards;
    else
        B = Breverse;
    end
    Load(i) = - B * thetadot + coulombicFriction(tau, thetadot);
end
end

function [xdot] = motorEQMotion(t, x)
% let the state be x = [theta; thetadot]
% then this function returns xdot = [thetadot; thetadotdot]

global J Bforwards Breverse

theta = x(1); thetadot = x(2);

if thetadot >= 0
    B = Bforwards;
else
   B = Breverse; 
end

T = motorTorque(motorCurrent(t));

%thetadotdot = T / J - B * thetadot / J + coulombicFriction(T, thetadot) / J;
thetadotdot = T / J;

xdot = [thetadot; thetadotdot];

end

function [It] = motorCurrent(t)
global Tobs Iobs
It = flatInterpolate(Tobs, Iobs, t);
end

function [torque] = motorTorque(Im)
% takes in motor current, Im
global Kt
torque = Kt * Im;
end

function [torque] = coulombicFriction(motorTorque, thetadot)
% coulombic friction cannont exceed the magnitude of motor torque and
% opposes attempted motion
% models friction as different in each direction
global Tfforwards Tfreverse StaticFrictionX

if (abs(thetadot) < 1e-8) % if not moving, friction opposes attempted motion
    % if motortorque does not exceed static friction
    if (motorTorque < Tfforwards && motorTorque >= 0) || (motorTorque > -Tfreverse && motorTorque < 0)
        FF = abs(motorTorque);
    elseif (motorTorque > Tfforwards)
        FF = Tfforwards;
    else
        FF = Tfreverse;
    end
    torque = - StaticFrictionX * FF * sign(motorTorque);
    
else % if we are moving oppose the actual motion with kinetic friciton
    if thetadot > 0
        torque = -Tfforwards;
    else
        torque = Tfreverse;
    end
end

end