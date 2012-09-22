function [phiMain phiLeak] = estimateFlux(I, x)
% given the current being supplied to the electromagnet and the
% distance x from the electromagnet to the top of the sphere,
% estmates the two components of the flux
% phiMain passes through the sphere while phiLeak does not

%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical constants
mu0 = 1.2566371e-6; % absolute permeability of vacum (H/m)
muair = 1.00000037 * mu0; % absolute permeability of air
% note that other permeabilities where already dropped out as being 
% identical to air
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical dimensions relavent to calculation

% electromagnet
Dem = 0.04602; % electromagnet diameter
Lem = 0.035; % electromagnet length
Aem = pi * (Dem / 2) * (Dem / 2); % effective area of electromagnet

% sphere
Ds = 0.0127; % sphere diameter
Ls = 0.0127; % length of path through sphere reflects assumption that 
% sphere is really a cylinder
As = pi * (Ds / 2) * (Ds / 2); % effective area of sphere taken as
% cross-section of cylinder contaiing sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approximations
N = 1000; % number of turns on the electromagnet TODO UNKNOWN

% point a is just outside the magnet on the upper edge
% point b is in the air next to the magnet at the level with the lowe edge
% point c is on the bottom of the sphere
% point d is on top of the sphere
% point e is on the bottom of the electromagnet
Lab = 1.2 * ((Dem / 2) + Lem); % path through air from a to b
Lbc = 1.2 * ((Dem / 2) + x + Ls); % path through air from b to c
Lbe = 1.2 * (Dem / 2);


% assumptions
% constant flux density
% conserved cross-sectional area - no fringing
% effective area from points a to b is equal to Aem
% effective area from point e to c is equal to As
% effective area from b to e is equal to Aem - As
% neglecting aluminum and copper in diagram because their permeabilities
%   are effectively equal to air
% Droping reluctance associated with steel sphere because musteel >> muair

%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute a few terms
C0 = (muair * (Aem - As) * Aem) / (Lbe * Aem + (Lem + Lab) * (Aem - As));
C1 = (Lem + Lab) / (muair * Aem);
C2 = Lbc / (muair * Aem);

% compute fluxes
phiMain = (N * I * (1 - C0 * C1)) / (C1 - C0 * C1 * C1 + C2 + x / (muair * As));
phiLeak = (N * I - phiMain * C1) * C0;

end

