function [] = parameterFittingTest()
% systme simulation, simulates the motor
% J*thetadotdot = sum(Forces) = torque - B * thetadot

global J B Kt Im Fk Fs

% note that these values are made up
J = 1; % moment of inertia of the rotor, constant
B = 0.1; % viscous damping 
Fk = 0.05; % torque exerted by coulombic friction, kinetic friction
Fs = 0.07; % torque exerted by static friction
Kt = 1;
Im = 0.5; % constant motor voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Controls

global CoulombicFrictionON
CoulombicFrictionON = 1;

ViscousDampingON = 1;
if (ViscousDampingON == 0)
    B = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let the state of the system 

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

xstart = [0; 0]; % initial angle, initial angular velocity
tend = 50;

[T, X] = ode45(@motorEQMotion, [0, tend], xstart, options);

I = motorTorque(motorCurrent(T));

Energy = 0.5 * J * X(:,2).^2;

%%%%%%%%%%%%%%%%
% Test parameter Estimation

for i = 1:length(T)
    if mod(i,10) == 0
        Ts(i/10) = T(i); Xs(i/10) = X(i,1); Is(i/10) = I(i);
    end
end

Pin = [J; B; Fk; Kt];
disp('giving points')
length(Ts)
Pc = parameterEstimation(0,Ts',Xs',Is',Pin);

disp('recieved P = ')
Pc

disp('actual P = ')
Pactual  = [J; B; Fk; Kt]

%%%%%%%%%%%%%%%%%%%%%%
% plot stuff
% figure;
% hold on;
% plot(T,X(:,1),'b')
% xlabel('Time (seconds)')
% ylabel('Angle (Radians)')
% title('Motor Angle vs Time')
% hold off;
% 
% figure
% hold on
% plot(T,X(:,2),'g')
% xlabel('Time (seconds)')
% ylabel('Angular Velocity (rads/s)')
% title('Motor Angular Velocity vs Time')
% hold off
% 
% figure 
% hold on;
% plot(T,I,'r')
% xlabel('Time (seconds)')
% ylabel('Motor Current (Amps)')
% title('Motor Current vs Time')
% hold off

% figure
% hold on;
% plot(T,X(:,1),'b')
% plot(T,X(:,2),'g')
% plot(T,I,'r')
% xlabel('Time')
% legend('Motor Position (radians)','Angular Velocity (rads/s)','Motor Current (Amps)')
% title('All')
% hold off
% 
% figure
% hold on;
% plot(T,Energy)
% xlabel('Time (seconds)')
% ylabel('Energy (Joules)')
% title('Energy vs Time')
% hold off

end

function [xdot] = motorEQMotion(t, x)
% let the state be x = [theta; thetadot]
% then this function returns xdot = [thetadot; thetadotdot]

global J B 
theta = x(1); thetadot = x(2);

T = motorTorque(motorCurrent(t));

thetadotdot = T / J - B * thetadot / J + coulombicFriction(T, thetadot);

xdot = [thetadot; thetadotdot];

end

function [current] = motorCurrent(t)
global Im
current = Im * sin(t);
end

function [torque] = motorTorque(Im)
% takes in motor current, Im
global Kt
torque = Kt * Im;
end

function [torque] = coulombicFriction(motorTorque, thetadot)
% coulombic friction cannont exceed the magnitude of motor torque and
% opposes attempted motion
global Fk Fs CoulombicFrictionON

if (CoulombicFrictionON == 0)
    torque = 0;
else
    if (abs(thetadot) < 1e-6) % if not moving, friction opposes attempted motion
        if abs(motorTorque) < Fs
            FF = abs(motorTorque); % and doesn't exceed the magnitude of the applied torque
        else
            FF = Fs;
        end
        torque = - FF * sign(motorTorque);
        
    else % if we are moving oppose the actual motion with kinetic friciton
        torque = - Fk * sign(thetadot);
    end
end

end
