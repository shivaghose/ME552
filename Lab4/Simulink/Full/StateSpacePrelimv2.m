% version 2 of the parameter initialization

% note that positive theta is counter clockwise when viewed from above
% positive alpha is clockwise when viewed from the end of the arm
% x axis points along pendulum
% z axis along arm towards motor
% y is consistent with right handed coordinate frame

% note that pendulum parameters have been obtained through a combination of
% solid works modeling to provide center of mass and moments of inertia
% and a rough eyeball fit of the friction and damping
% arm damping and friction terms are carried over from the previous
% experiment but this may not be quite the case with an off balance torque
% on the shaft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical parameters

% Motor and Driver Parameters
Kt = 0.0314499;     % motor torque constant
Kb = 0.04;          % yes technically these should be equal, this is based on 600 rad/s 
                    % for a supply voltage of 24 volts only effects driver /
                    % motor model
Rm = 1.85;          % motor resistance based on data sheet
Vsupply = 24;       % supp,y voltage
Ka = 1;             % nominal driver voltage to current ratio


% lengths
l1 = 0.10825640;
l2c = 0.0856786;

% mass
m2 = 0.0820142;

% moments of inertia
Ix2 = 8.1e-6;
Iy2 = 3.711e-4;
Iz2 = 3.646e-4;
Ixz2 = 2.03e-5;
I1zt = 6.78e-4 + 3.50514e-5;    % arm moment of inertia obtained from solid works
                                % and rotor inertia obtained from earlier
                                % experiments

% damping 
ba = 6.5E-6;        % alpha / pendulum    
btf = 1.08586e-5;   % forwards theta / arm  
btr = 4.8593e-5;    % reverse theta

% friction
Tca = 7e-6;             % alpha / pendulum
Tctf = 0.022439;        % forwards theta / arm
Tctr = 0.0143096;       % reverse theta
staticThreshold = 1e-6; % below this velocity we aren't moving
staticRatio = 1;        % multiplier for static friction over normal friction
% alternate friction quantities which do not illustrate a limit cycle
%Tctf = 0.0001; Tctr = 0.0001;

% gravity
g = 9.81; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transfer functions
% this linear form requires us to neglect friction and assume isotropic
% damping

% constants
C1 = Iz2+m2*l2c*l2c;  C2 = Iy2+m2*l2c*l2c;  C3 = m2*l1*l2c-Ixz2;
C4 = l2c*m2*g;        C5 = I1zt+Ix2+m2*l1*l1;   C6 = C5*C1-C3*C3;

btavg = (btf + btr)/2; 
bt = btavg;

% with Damping
DenB = [(C5*C1-C3*C3), (C1*btavg+C5*ba), (ba*btavg-C4*C5), (-C4*btavg), 0];
NumBT = [C1, ba, -C4];
NumBA = [-C3, 0, 0]; 

tfTB = tf(NumBT,DenB);
tfAB = tf(NumBA,DenB);

% without damping
Den = [(C5*C1-C3*C3), 0, -C4*C5, 0, 0];
NumT = [C1, 0, -C4];
NumA = [-C3, 0, 0];

tfT = tf(NumT,Den);
tfA = tf(NumA,Den);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Space
% where x = [theta; thetadot; alpha; alphadot];
% xdot = A*x + B*T

A = [  0,  1,              0,              0;
        0,  (-C1*btavg/C6), (-C3*C4/C6),    (C3*ba/C6);
        0,  0,              0,              1;
        0,  (C3*btavg/C6),  (C4*C5/C6),     (-C5*ba/C6)];
    
B = [0; C1/C6; 0; -C3/C6];    

t = 0:.01:9.99;
C = [0 0 1 0];
D=0;

sys = ss(A,B,C,D);
x0 = [0 0 .1 0]';
u=zeros(size(t));

% 10 10 100 40 R = 5
Q = diag([10 40 120 60]);
R = 6;

K = lqr(A,B,Q,R)
K1 = K(1); K2 = K(2); K3 = K(3); K4 = K(4);

Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;

sys_cl = ss(Ac,Bc,Cc,Dc);

% [y,x] = lsim(A,B,C,D, u, t, x0);
% plot(t,x)
% title('Uncontrolled');
% figure
% [Y,T,X] = lsim(sys_cl,u,t, x0);
% plot(T,Y);
% title('Controlled')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Controller Parameters
CONTROLLER_SWITCH_POINT = .30;  % switch over to balancing controller here
tau = .01;                      % time constant for derivative filter
    
    