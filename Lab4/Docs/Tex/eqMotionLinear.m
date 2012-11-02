function [tdot_, Tdotdot, adot_, Adotdot] = eqMotionLinear(Tm, theta,tdot,alpha,adot, m2, I1zt, Ix2, Iy2, Iz2, Ixz2, g, l1, l2c, ba, bt)
tdot_ = tdot;
adot_ = adot;

%%%%%%%%%%%%
% Linear Damping
% Neglect Friction

C1 = Iz2+m2*l2c*l2c;  C2 = Iy2+m2*l2c*l2c;   C3 = m2*l1*l2c-Ixz2;  C4 = l2c*m2*g;

C5 = I1zt + Ix2 + m2*l1*l1;  C6 = C5*C1-C3*C3;

Tdotdot = (1/C6)*(C3*ba*adot - C3*C4*alpha - C1*bt*tdot + C1*Tm);

Adotdot = (1/C6)*(-ba*adot*C5 + bt*tdot*C3 + alpha*C4*C5 - C3*Tm);
end