function [stateC] = enchancedParameterEstimation(Tin,Thetain,Iin, state0, lockedPositions)
%Oh god the state so many terms
% where state0 = [J; Kt; Bf; Br; Tff; Tfr]
% Bf and Tff are the viscous and coulomb damping forces applied when
% rotating or attempting to rotate in the positive direction
% Br and Tfr for the reverse direction
% lockedPositions is a vector that contains non zero elements in positions
% that will not be optimized

global state LPos numUPos staticRatio staticThreshold Tobs Iobs Thetaobs numPoints
assert((length(state) == 6),'Invalid state length');
assert((length(state) == length(lockedPositions)),'should have been the same length');

% ode45 sort of works, but is very slow, static friction really gives it a
% hard time

numPoints = 2*10^4;
state = state0;
LPos = lockedPositions;

% how many state variables are we actually going to regress on
numUPos = length(state);
for i = 1:length(LPos)
    if LPos(i)
       numUPos = numUPos - 1; 
    end
end

staticRatio = 1.0;  % static friction assumed to be 10% larger than sliding friction to avoid stuttering problem
staticThreshold = 1e-6; % below this we aren't moving
Tobs = Tin; Thetaobs = Thetain; Iobs = Iin;

xstart = [0; 0];
stateC = state0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial plot to make sure the shooting works correctly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T, X] = shootMotor(xstart,stateC);
figure
hold on;
plot(T,X(:,1),'b');
plot(Tobs,Thetaobs,'k');
title('Initial Estimate curve');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
hold off;

dPtol = 1e-6;
dP = inf;
iter = 0;
maxiter = 1000;
lastE = inf;
lastS = stateC;

while ((max(abs(dP)) >= dPtol) && (iter < maxiter))
    disp(iter)
    E = computeError(xstart,stateC);
   
    % call it quits if the error increased, w
    magE = norm(E)
    if magE > lastE
        stateC = lastS;
        break
    end
    lastE = magE;
    lastS = stateC;
    
    J = jacobianShootMotor(xstart,stateC);
    
    ds = J \ E;
    
    % going 100% of the dP is causing trouble, going less distance should
    % help
    stateC = updateState(stateC,0.1 .* ds);
    
    iter = iter + 1;
end

E = computeError(xstart,stateC);

disp('magnitude of final error =')
norm(E)

%%%%%%%%%%%%%%%%%%
% Make some plots
%%%%%%%%%%%%%%%%%%
[T, X] = shootMotor(xstart,stateC);

figure
hold on;
plot(T,X(:,1),'b');
plot(Tobs,Thetaobs,'k');
title('Fit curve');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
hold off;

end

function [J] = jacobianShootMotor(xstart,state0)
% compute the jacobian of the error vector in the neighborhood of the
% values P0

global Tobs numUPos LPos

[T0,X0] = shootMotor(xstart,state0);

% matrix for storing perturbed values
pX1 = zeros(length(Tobs),numUPos);
deltaP = zeros(numUPos,1);

dP = 1e-2;

p = 1;
for s = 1:length(state0)
    
    if ~LPos(s)
        stateP = state0;
        deltaP(p) = dP * abs(state0(s));
        %deltaP(p) = 1e-6;
        stateP(s) = stateP(s) + deltaP(p) ;
        [pT, pX] = shootMotor(xstart,stateP);
        pX1(:,p) = pX(:,1); % we only care about storing the positions
        % don't discretize yet, because that would mess up the jacobian
        p = p + 1;
    end
end

J = zeros(length(T0),numUPos);
for i = 1:length(T0)
   for j = 1:numUPos
       J(i,j) = (pX1(i,j) - X0(i)) / deltaP(j);
   end
end

end

function [E] = computeError(xstart,state)
global Thetaobs

[T, X] = shootMotor(xstart,state); % don't need to interpolate, now runs exactly at the sampled points

E = Thetaobs - X(:,1); 
end

function [T, X] = shootMotor(xstart,state_)

global state Tobs numPoints

state = state_;

tstart = Tobs(1);
tend = Tobs(end);

%options = odeset('RelTol',1e-7,'AbsTol',1e-7);
%[T, X] = ode45(@motorEQMotion, Tobs,xstart,options);

T = linspace(tstart,tend, numPoints);
[X] = ode5(@motorEQMotion, T, xstart);
Xint = interpolate(T,X,Tobs);

X = Xint;
T = Tobs;
end

function [snew] = updateState(s, ds)
% can't just do s = s + ds because the actual state and the number of
% parameters we are regressing aren't equal

global LPos

snew = s;
j = 1;
for i = 1:length(s) % this is the full length state
    if ~LPos(i) % if not locked then can update
        snew(i) = snew(i) + ds(j);
        j = j + 1;
    end
end
end

function [xdot] = motorEQMotion(t, x)
% x = [theta; thetadot]

thetadot = x(2); % theta = x(1);

global state staticThreshold
J = state(1);

% need to choose the right value for B based on which direction we are
% turning in
if (thetadot > staticThreshold)
    B = state(3); % Bf = state(3);
elseif thetadot < - staticThreshold
    B = state(4); % Br = state(4);
else
    B = 0; % to avoid spurious forces when we aren't moving
end

T = motorTorque(motorCurrent(t));

thetadotdot = T / J - B * thetadot / J + coulombFriction(T, thetadot) / J;

xdot = [thetadot; thetadotdot];
end

function [It] = motorCurrent(t)
global Tobs Iobs
It = flatInterpolate(Tobs,Iobs,t);
end

function [torque] = motorTorque(Im)
global state  %Kt = state(2);
torque = state(2) * Im;
end

function [torque] = coulombFriction(motorTorque, thetadot)
% coulombic friction cannont exceed the magnitude of motor torque and
% opposes attempted motion
% models friction as different in each direction
% forces are not equal in each direction
global state staticThreshold staticRatio
Tff = state(5); Tfr = state(6);

if (abs(thetadot) < staticThreshold) % if not moving, friction opposes attempted motion
    % if motortorque does not exceed static friction
    if ((motorTorque < Tff) && (motorTorque >= 0)) || ((motorTorque > -Tfr) && (motorTorque < 0))
        FF = abs(motorTorque);
    elseif (motorTorque > Tff)
        FF = Tff;
    else
        FF = Tfr;
    end
    torque = - staticRatio * FF * sign(motorTorque);
    
else % if we are moving oppose the actual motion with kinetic friciton
    if thetadot > staticRatio
        torque = -Tff;
    else
        torque = Tfr;
    end
end

end

function [Xout] = interpolate(Tin,Xin,Tout)
% linear interpolation of the desired Tout values

Xout = zeros(length(Tout),1);

i = 1; % keep track of positon in the original
for o = 1:length(Tout)
    while (i <= length(Tin)) && (Tout(o) > Tin(i))
        i = i + 1;
    end
    if (i == 1) % if we asked for a T before the first T
        dt = Tin(2) - Tin(1);
        slope = (Xin(2) - Xin(1)) / dt;
        dt = Tout(o) - Tin(1);
        Xout(o) = Xin(1) + slope * dt;
    elseif (i > length(Tin))
        % linearly extrapolate - subject to noise, be careful
        dt = Tin(end) - Tin(end - 1);
        slope = (Xin(end) - Xin(end - 1)) / dt;
        dt = Tout(o) - Tin(end);
        Xout(o) = Xin(end) + slope * dt;
    else
        dt = Tin(i) - Tin(i - 1);
        slope = (Xin(i) - Xin(i - 1)) / dt;
        dt = Tout(o) - Tin(i - 1);
        Xout(o) = Xin(i - 1) + slope * dt;
    end
end

end