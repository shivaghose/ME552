function [PBest] = parameterSweep(Tin, Thetain, Iin, Pin)
% sweep parameters in a wide range to reduce error

global Tobs Thetaobs Iobs
Tobs = Tin; Thetaobs = Thetain; Iobs = Iin; P0 = Pin;
xstart = [0; 0];

deltaP = 0.4 .* P0; % use 0.1% of initial value for the sweep sort of like resolution
Range = 2;

% try the initial
PBest = P0;
EBest = norm(computeError(xstart,PBest));

% have 4 parameters to sweep
iter = 0;
for a = -Range : Range
   for b = -Range : Range
      for c = -Range : Range 
          for d = -Range : Range
              iter = iter + 1
              
              Ptest = P0;
              Ptest(1) = Ptest(1) + a * deltaP(1);
              Ptest(2) = Ptest(2) + b * deltaP(2);
              Ptest(3) = Ptest(3) + c * deltaP(3);
              Ptest(4) = Ptest(4) + d * deltaP(4);
              
              if ((Ptest(1) < 0) || (Ptest(2) < 0) || (Ptest(3) < 0) || (Ptest(4) < 0))
                  continue;
              end
              
              Etest = norm(computeError(xstart,Ptest));
              % if it got better keep it
              if (Etest < EBest)
                  EBest = Etest;
                  PBest = Ptest;
              end
          end
      end
   end
end

[T, X] = shootMotor(xstart,PBest);

figure
hold on;
plot(T,X(:,1),'b');
plot(Tobs,Thetaobs,'ko');
title('Fit curve');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
hold off;

end

% computes errorvector between our curve and the observed values
function [E] = computeError(xstart,P)

global Tobs Thetaobs

[T,X] = shootMotorKnownT(xstart,P,Tobs);

E = Thetaobs - X(:,1);
end

function [J] = jacobianShootMotor(xstart,P0)
% compute the jacobian of the error vector in the neighborhood of the
% values P0

global P Tobs

[T0,X0] = shootMotorKnownT(xstart,P0,Tobs);

% matrix for storing perturbed values
pX1 = zeros(length(T0),length(P0));

deltaP = 1e-9;

for p = 1:length(P0)
    pP = P0;
    pP(p) = pP(p) + deltaP;
    [pT, pX] = shootMotorKnownT(xstart,pP,Tobs);
    pX1(:,p) = pX(:,1); % we only care about storing the positions
    % don't discretize yet, because that would mess up the jacobian
end

J = zeros(length(T0),length(P0));
for i = 1:length(T0)
   for j = 1:length(P0)
       J(i,j) = (pX1(i,j) - X0(i)) / deltaP;
   end
end

end

% where x0 are the initial conditions initial angle, initial angular
% velocity, and P_ is the value of motor parameters to use
function [T, X] = shootMotor(xstart,P_)

global Tobs P

P = P_;

% determine starting and ending times from our data
tstart = 0;
tend = Tobs(end) + 0.25; % go quarter of a second longer past the end

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

[T, X] = ode45(@motorEQMotion, [tstart, tend], xstart, options);

end

function [T, X] = shootMotorKnownT(xstart,P_,Tin)
% if we want our solution values at the Tin times

global P

P = P_;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

[T, X] = ode45(@motorEQMotion, Tin, xstart, options);

end

function [xdot] = motorEQMotion(t, x)
% let the state be x = [theta; thetadot]
% then this function returns xdot = [thetadot; thetadotdot]

global P
J = P(1); B = P(2);
theta = x(1); thetadot = x(2);

T = motorTorque(motorCurrent(t));

thetadotdot = T / J - B * thetadot / J + coulombicFriction(T, thetadot);

xdot = [thetadot; thetadotdot];

end

function [It] = motorCurrent(t)
global Tobs Iobs
It = flatInterpolate(Tobs, Iobs, t);
end

function [torque] = motorTorque(Im)
% takes in motor current, Im
global P
Kt = P(4);
torque = Kt * Im;
end

function [torque] = coulombicFriction(motorTorque, thetadot)
% coulombic friction cannont exceed the magnitude of motor torque and
% opposes attempted motion
global P

Ff = P(3);

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