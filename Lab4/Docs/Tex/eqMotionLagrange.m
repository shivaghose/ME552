function [tdot_, Tdotdot, adot_, Adotdot]  = eqMotionLagrange(motorTorque, frictionOn , dampingOn ,theta, tdot, alpha, adot, Iz2,l2c,Iy2,l1,Ixz2,I1zt,Ix2,m2,g,ba,btf,btr,Tca,Tctf,Tctr)
% we currently assume differential coulomb friction and damping for theta,
% but not for alpha
%global m2 I1zt Ix2 Iy2 Iz2 Ixz2 g l1 l2c ba btf btr

% Iz2 = 3.646e-4;
% l2c = 0.0856786;
% Iy2 = 3.711e-4;
% l1 = 0.1016;
% Ixz2 = 2.03e-5;
% I1zt = 6.78e-4;
% Ix2 = 8.1e-6;
% m2 = 0.0820142;
% g = 9.81; ba = 6.5E-6; btf = 1.08586e-5; btr = 4.8593e-5;
% %Tca = 7e-6; Tctf = 0.022439; Tctr = 0.0143096; 
% Tca = 7e-6; Tctf = 0.00; Tctr = 0.00; 

% theta = x(1); tdot = x(2); ct = cos(theta); st = sin(theta);
% alpha = x(3); adot = x(4); ca = cos(alpha); sa = sin(alpha);
tdot_ = tdot;
adot_ = adot;
%ct = cos(theta); st = sin(theta);
ca = cos(alpha); sa = sin(alpha);
%%%%%%%%%%%
% Motor
Tm = motorTorque;  % placeholder for motor power term

%%%%%%%%%%%%
% Damping
% model damping differently forwards and backwards for theta
if dampingOn
    if (tdot > 0)
        Bt = btf;  % baf opposes forwards motion for theta
    elseif (tdot < 0)
        Bt = btr; % bar opposes reverse motion for theta
    else
        Bt = 0;
    end
    Ba = ba;
else
    Bt = 0;
    Ba = 0;
end
% for alpha only a single viscous damping term

%%%%%%%%%%%%%%
% Friction
if frictionOn
    TCt = coulombFriction(Tctf, Tctr, Tm, tdot);
    Tg = m2*g*l2c*sa; % torque exerted by gravity on pendulum
    TCa = coulombFriction(Tca, Tca, Tg, adot);
else 
    TCa = 0;
    TCt = 0;
end


C1 = Iz2+m2*l2c*l2c;  C2 = Iy2+m2*l2c*l2c;   C3 = m2*l1*l2c-Ixz2;  C4 = l2c*m2*g;

D = (I1zt + sa*sa*C2 + Ix2*ca*ca + m2*l1*l1)*C1 - ca*ca*C3*C3;

Tdotdot = (1/D)*(-2*tdot*adot*sa*ca*(C2-Ix2)*C1 + adot*adot*sa*C3*C1 ...
    - tdot*tdot*sa*ca*ca*(C2-Ix2)*C3 + adot*Ba*ca*C3 - sa*ca*C3*C4 - TCa*ca*C3 ...
    - Bt*tdot*C1 + Tm*C1 + TCt*C1);

Adotdot = (1/C1)*(-Tdotdot*ca*C3 + tdot*tdot*sa*ca*(C2-Ix2) + C4*sa - Ba*adot + TCa);

end
% xdot = [tdot; Tdotdot; adot; Adotdot];

%-------------------------------------------------
function [torque] = coulombFriction(Tff, Tfr, externalTorque, thetadot)
% this friction model is making an assumption that isn't strictly true, we
% have static friction, but only the external torque exerted on each link
% will be considered.  This is not strictly correct, really we should
% consider the constraint forces as well
% it might not hold with the strength that it should

staticThreshold = 1e-6;
staticRatio = 1;

if (abs(thetadot) < staticThreshold) % if not moving, friction opposes attempted motion
    % if external torque does not exceed static friction
    if ((externalTorque < Tff) && (externalTorque >= 0)) || ((externalTorque > -Tfr) && (externalTorque < 0))
        FF = abs(externalTorque);
    elseif (externalTorque > Tff)
        FF = Tff;
    else
        FF = Tfr;
    end
    torque = - staticRatio * FF * sign(externalTorque);
    
else % if we are moving oppose the actual motion with kinetic friciton
    if thetadot > staticThreshold
        torque = -Tff;
    else
        torque = Tfr;
    end
end
end
%------------------------------------------------
