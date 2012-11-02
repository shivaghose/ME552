function [tdot_, Tdotdot, adot_, Adotdot] = eqMotionLinearDown(Tm, theta,tdot,alpha,adot, m2, I1zt, Ix2, Iy2, Iz2, Ixz2, g, l1, l2c, ba, bt)
% linearized about alpha = pi, theta = 0
tdot_ = tdot;
adot_ = adot;

%%%%%%%%%%%%
% Linear Damping
% Neglect Friction

C1 = Iz2+m2*l2c*l2c;  C2 = Iy2+m2*l2c*l2c;   C3 = m2*l1*l2c-Ixz2;  C4 = l2c*m2*g;

C7 = (I1zt + C2*(pi*pi - 2*pi*alpha) + Ix2 + m2*l1*l1)*C1 - C3*C3;

Tdotdot = (1/C7)*(C3*C4*(alpha - pi) - C3*ba*adot - C1*bt*tdot + C1*Tm);

Adotdot = (1/C1)*(Tdotdot*C3 + C4*pi - C4*alpha - ba*adot);
end