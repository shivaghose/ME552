function [energy, KE, PE]  = totalEnergy(theta, tdot, alpha, adot, Iz2,l2c,Iy2,l1,Ixz2,I1zt,Ix2,m2,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute potential energy of the system only cares about the pendulum
ca = cos(alpha); sa = sin(alpha);
ct = cos(theta); st = sin(theta);

PE = (ca - 1)*l2c*m2*g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the kinetic energy of the the system
thetadot = tdot;  alphadot = adot;

% rotational component for arm
KE1 = 0.5*I1zt*thetadot*thetadot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translational component for pendulum
J = zeros(3,2);

J(1,1) = -l2c*sa*ct - l1*st;
J(1,2) = -l2c*st*ca;

J(2,1) = -l2c*sa*st + l1*ct;
J(2,2) = l2c*ct*ca;

J(3,2) = -l2c*sa;

% now that we have the jacobian
V = J*[thetadot; alphadot];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KE2T = 0.5*m2*(V'*V);

% rotaional component
KE2R = 0.5*(thetadot*thetadot*(Ix2*ca*ca + Iy2*sa*sa) - 2*Ixz2*thetadot*alphadot*ca ...
    + Iz2*alphadot*alphadot);

KE = KE1 + KE2T + KE2R;

energy = KE + PE;