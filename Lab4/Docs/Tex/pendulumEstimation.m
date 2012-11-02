function [stateC] = pendulumEstimation(Tin, Alphain, xstart, state0, lockedVariables)
% only consider the pendulum by itself, because it is much simpler

% state for simulation
% x = [alpha; alphadot]

% state for optimization
% state = [Iz2 m2 l2c ba Tca]

assert((length(state0) == 5), 'Invalid state length, expected 5');
assert((length(state0) == length(lockedVariables)),'should have been the same length');
assert((length(xstart) == 2), 'We need initial conditions');

global state g staticRatio staticThreshold windowSize TIMES
state = state0;
g = 9.81;
staticRatio = 1.1; staticThreshold = 1e-8; windowSize = 20; TIMES = 3;

global Tobs Alphaobs LPos numUPos 
Tobs = Tin;
Alphaobs = Alphain;
LPos = lockedVariables;

% how many state variables are we actually going to regress on
numUPos = length(state);
for i = 1:length(LPos)
    if LPos(i)
       numUPos = numUPos - 1; 
    end
end

% nominal state
%state = [3.646e-4; 0.0820142; 0.0856786; 1e-5; 1e-5]; 
%xstart = [pi/2; 0];

%%%%%%%%%%%%%%%
% Initial
%%%%%%%%%%%%%%%
[T, X] = shootCompare(xstart, state);
stateC = state0;

plotEnvelopes(T,X)
plotCompare(T,X);
%plotEnergy(T,X);

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
    
    J = jacobianShootPendulum(xstart,stateC);
    
    ds = J \ E;
    
    % going 100% of the dP is causing trouble, going less distance should
    % help
    stateC = updateState(stateC,ds);
    
    iter = iter + 1;
end

E = computeError(xstart,stateC);

disp('magnitude of final error =')
norm(E)

%%%%%%%%%%%%%
% Final 
[T, X] = shootCompare(xstart, state);

plotComparey(T,X);
plotEnergy(T,X);

end

function [J] = jacobianShootPendulum(xstart,state0)
% compute the jacobian of the error vector in the neighborhood of the
% values P0

global Tobs numUPos LPos

[T0,X0] = shootCompare(xstart,state0);

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
        [pT, pX] = shootCompare(xstart,stateP);
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
global Alphaobs

[T, X] = shootCompare(xstart,state); % don't need to interpolate, now runs exactly at the sampled points

E = Alphaobs - X(:,1); 
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


function [T, X] = shootNormal(xstart,tend, state_)
global state
state = state_;

options = odeset('RelTol',1e-6,'AbsTol',1e-6);
[T, X] = ode45(@EQMotion,[0 tend],xstart,options);
end

function [T, X] = shootCompare(xstart, state_)
global Tobs state TIMES
state = state_;
%options = odeset('RelTol',1e-6,'AbsTol',1e-6);
%[T, X] = ode45(@EQMotion,Tobs,xstart,options);
Tint = interpolateTime(Tobs,TIMES);
Xint = ode5(@EQMotion,Tint,xstart);
X = reduceTime(Xint,TIMES);
T = Tobs;
end

function [Tout] = interpolateTime(T,times)
Tout = zeros((length(T)*times - (times-1)), 1);

for i = 1:(length(T)-1)
    dT = T(i+1) - T(i); 
   for j = 0:(times-1)
       Tout((i*times-(times-1))+j) = T(i) + j*dT/times;
   end
end

Tout(end) = T(end);

end

function [Tout] = reduceTime(T,times)
% times is the times that was applied with interpolate
Tout = zeros((length(T)+(times-1))/times,1);

for i = 1:length(Tout)
   Tout(i) = T(i*times-(times-1)); 
end
end

function [xdot] = EQMotion(t,x)
% equation of motion for just the pendulum not all that many parameters
global state g

Iz2 = state(1); m2 = state(2); l2c = state(3); ba = state(4); Tca = state(5);

alpha = x(1); alphadot = x(2); sa = sin(alpha); 

Tg = m2*g*l2c*sa;
Tc = coulombFriction(Tca, Tg, alphadot);

D = Iz2 + m2*l2c*l2c;

alphadotdot = (1/D)*(-ba*alphadot + m2*g*l2c*sa + Tc);

xdot = [alphadot; alphadotdot];
end


function [torque] = coulombFriction(Tf, externalTorque, thetadot)
% this friction model is making an assumption that isn't strictly true, we
% have static friction, but only the external torque exerted on each link
% will be considered.  This is not strictly correct, really we should
% consider the constraint forces as well
% it might not hold with the strength that it should
global staticThreshold staticRatio

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

function [Xout] = interpolate(Tin,Xin,Tout)
% linear interpolation of the desired Tout values

Xout = zeros(length(Tout),size(Xin,2));

i = 1; % keep track of positon in the original
for o = 1:length(Tout)
    for j = 1:size(Xin,2)
        while (i <= length(Tin)) && (Tout(o) > Tin(i))
            i = i + 1;
        end
        if (i == 1) % if we asked for a T before the first T
            dt = Tin(2) - Tin(1);
            slope = (Xin(2,j) - Xin(1,j)) / dt;
            dt = Tout(o) - Tin(1);
            Xout(o,j) = Xin(1,j) + slope * dt;
        elseif (i > length(Tin))
            % linearly extrapolate - subject to noise, be careful
            dt = Tin(end) - Tin(end - 1);
            slope = (Xin(end,j) - Xin(end-1,j)) / dt;
            dt = Tout(o) - Tin(end);
            Xout(o,j) = Xin(end,j) + slope * dt;
        else
            dt = Tin(i) - Tin(i - 1);
            slope = (Xin(i,j) - Xin(i-1,j)) / dt;
            dt = Tout(o) - Tin(i-1);
            Xout(o,j) = Xin(i-1,j) + slope * dt;
        end
    end
end
end


function [] = plotCompare(T, X)
global Tobs Alphaobs
% plot the system
figure
hold on
plot(T,X(:,1),'k');
plot(Tobs,Alphaobs,'b')
legend('Alpha','Alpha Observed');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
title('Passive Pendulum Behavior')
end

function [] = plotTrajectory(T, X)
% plot the system
figure
hold on
plot(T,X(:,1),'k');
plot(T,X(:,2),'b');
legend('Alpha','Alphadot');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
title('Passive Pendulum Behavior')
end

function [] = plotEnergy(T,X)
% plots energy of the system

PE = potentialEnergy(X);
KE = kineticEnergy(X);
Total = PE + KE;

figure;
hold on;
plot(T,PE,'g');
plot(T,KE,'r');
plot(T,Total,'k');
xlabel('Time (seconds)');
ylabel('Energy');
legend('Potential','Kinetic','Total');
title('Energy');
hold off

end

function [] = plotEnvelopes(T,X)
global Tobs Alphaobs windowSize
[Tout,UA,LA] = extractEnvelop(Tobs,Alphaobs,windowSize);
[Tout,UX,LX] = extractEnvelop(T,X,windowSize);
figure
hold on
plot(Tout,UA,'b')
plot(Tout,LA,'b')
plot(Tout,UX,'k')
plot(Tout,LX,'k')
xlabel('Time (seconds)');
ylabel('Angle (radians');
title('Trajectory Envelope');
hold off
end

function KE = kineticEnergy(X)
% what is the kinetic energy of the system
global state
Iz2 = state(1); m2 = state(2); l2c = state(3);

KE = 0.5 * (Iz2 + m2*l2c*l2c) * X(:,2) .* (X(:,2));

end

function PE = potentialEnergy(X)
% what is the potential energy of the system
global state g
m2 = state(2); l2c = state(3);

PE = l2c*m2*g*(cos(X(:,1)) - 1);

end

function [] = test(T)
last = -1;
for i = 2:length(T)
    if (T(i) - T(i - 1)) < 0
        disp('we went down!')
        i
    end
end
end