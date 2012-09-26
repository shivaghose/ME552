function [dRdd] = radianceDerivative(rbeam, rsphere, d)
% compute derivative of the radiance curve at a particular distance

assert((rbeam >= 0) && (rsphere >= 0), 'Expected non-negative numbers')

if ((rsphere + abs(d)) <= rbeam)
    % handle the block is completely in the beam and the beam is bigger
    dRdd = 0;
    
elseif ((rbeam + abs(d)) <= rsphere)
    % handle case where block is larger or same size than beam and completely in the way
    dRdd = 0;
    
elseif (abs(d) < (rbeam + rsphere))
    Abeam = pi * rbeam * rbeam;
    
    % y stuff
    y = (d*d + rbeam*rbeam - rsphere*rsphere) / (2 * d);
    dydd = 1/2 - (rbeam*rbeam)/(2*d*d) + (rsphere * rsphere) / (2*d*d)
    
    % intersection distances
    dbeam = y;
    ddbeam = dydd
    dsphere = d - y;
    ddsphere = 1 - dydd
    
    % occluded areas
    AObeam = cSarea(rbeam,dbeam)
    dAObeam = cSDarea(rbeam, dbeam, ddbeam)
    AOsphere = cSarea(rsphere,dsphere)
    dAOsphere = cSDarea(rsphere, dsphere, ddsphere)
    
    % total overlap
    AOverlap = AObeam + AOsphere;
    dAOverlap = dAObeam + dAOsphere;
    
    % change in radiance fraction
    dRdd = dAOverlap / Abeam;
else
    % handle if they aren't even touching
    dRdd = 0;
end

end

function [area] = cSarea(r, h) 
% computes the area of the circular segment where the upper boundary is a
% circular arc and the lower boundary is a cord of the circle with radius r
% h is the perpendicular distance from the center to the cord

area = r * r * acos(h / r) - h * sqrt(r * r - h * h);

end

function [dA] = cSDarea(r, h, dh)
dA = (- (r / sqrt(1 - (h/r)^2)) - sqrt(r*r - h*h) + ((h*h) / sqrt(r*r - h*h)))*dh;
end