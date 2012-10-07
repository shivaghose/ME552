function [Pc] = parameterEstimation(filename,Tin,Thetain,Iin,Pin)
% Let P be the vector of model parameter values
% P = [J; B; Ff; Kt]
% J: rotor moment of inertia
% B: viscous damping
% Fs: force of friction, not modeling static friction seperately
% Kt: motor torque constant

global Tobs Thetaobs Iobs

if filename ~= 0
    % read data from file
    VAL = xlsread(filename);
    Tobs = VAL(:,1); Iobs = VAL(:,2); Thetaobs = VAL(:,3);
    P0 = [8.5e-6; 3.7e-6; (5.6e-3)/2; 4.24e-2]; % estimated values from the datasheet
else
    Tobs = Tin; Thetaobs = Thetain; Iobs = Iin; P0 = Pin;
end
xstart = [0; 0];
Pc = P0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial plot to make sure the shooting works correctly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[T, X] = shootMotor(xstart,Pc);
figure
hold on;
plot(T,X(:,1),'b');
plot(Tobs,Thetaobs,'ko');
title('Initial Estimate curve');
xlabel('Time (seconds)')
ylabel('Angle (radians)')
hold off;

dPtol = 1e-6;
dP = inf;
iter = 0;
maxiter = 1000;
lastE = inf;
lastP = Pc;

while ((max(abs(dP)) >= dPtol) && (iter < maxiter))
    disp(iter)
    E = computeError(xstart,Pc);
   
    % call it quits if the error increased, w
    magE = norm(E)
    if magE > lastE
        Pc = lastP;
        break
    end
    lastE = magE;
    lastP = Pc;
    
    J = jacobianShootMotor(xstart,Pc);
    
    dP = J \ E;
    
    % going 100% of the dP is causing trouble, going less distance should
    % help
    Pc = Pc + 0.1*dP;
    
    iter = iter + 1;
end

E = computeError(xstart,Pc);

disp('magnitude of final error =')
norm(E);

%%%%%%%%%%%%%%%%%%
% Make some plots
%%%%%%%%%%%%%%%%%%
[T, X] = shootMotor(xstart,Pc);

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

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[T, X] = ode45(@motorEQMotion, [tstart, tend], xstart, options);

end

function [T, X] = shootMotorKnownT(xstart,P_,Tin)
% if we want our solution values at the Tin times

global P

P = P_;

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

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

if (abs(thetadot) < 1e-8) % if not moving, friction opposes attempted motion
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