function [irradianceRatio] = irradianceEstimate(rbeam, rblock, d)
% computes the fraction of irradiance as a function of the radii of the two
% intersecting circles and the distance between their centers, d
% This can be used to model the irradiance of our phototransistor if we
% make several assumptions:
% - ignore diffraction
% - assume parallel light rays from LED to photo transistor so that we can
% model the * note we should use the minimum radius of the LED and the
% phototransistor - since - the parallel ray assumption means that distance
% from the ball to the LED doesn't matter

% let r1 be the beam radius and r2 be the blocking radius

% details at http://mathworld.wolfram.com/Circle-CircleIntersection.html

% estimated rblock = 6.35 mm
% estimated rbeam as = 2.3650 mm radius of the phototransistor

% regressed rbeam = 1.5048

assert((rbeam >= 0) && (rblock >= 0), 'The radii should be non-negative');

r1 = rbeam;
r2 = rblock;

if ((r1 > r2) && ((r2 + abs(d)) < r1))
    % handle the block is completely in the beam and the beam is bigger
    A1 = pi * r1 * r1;
    A2 = pi * r2 * r2;
    
    irradianceRatio = (A1 - A2) / A1;
    
elseif ((r2 >= r1) && ((r1 + abs(d)) <= r2))
    % handle case where block is larger or same size than beam and completely in the way
    irradianceRatio = 0;
        
elseif (abs(d) < (r1 + r2))
    % handle case where two circles intersect each other
    x = (d * d + r1 * r1 - r2 * r2) / (2 * d);
    
    d1 = x;
    d2 = (d * d - r1 * r1 + r2 * r2) / (2 * d);
    
    A = cSarea(r1, d1) + cSarea(r2, d2);
    
    Amax = pi * r1 * r1;
    
    irradianceRatio = (Amax - A) / Amax;

else 
    % handle case with no overlap
    irradianceRatio = 1;
end

end

function [area] = cSarea(r, h) 
% computes the area of the circular segment where the upper boundary is a
% circular arc and the lower boundary is a cord of the circle with radius r
% h is the perpendicular distance from the center to the cord

area = r * r * acos(h / r) - h * sqrt(r * r - h * h);

end