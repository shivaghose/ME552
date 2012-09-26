function [estRbeam estRblock] = estimateRbeamRblock()
% a bit of non-linear least squares to choose rbeam that best explains the
% data that we've seen

% initial estimate for R = [rbeam, rblock];
r = [2.3650; 6.35];

dRtol = 1e-9;
maxiter = 1000;

dR = Inf;
iter = 0;

% our observed data
displacement = [3.49; 4.45; 4.927; 5.35; 5.4; 5.85; 6.35; 6.85; 7.3; 7.78; 9.21];
VoltageFraction = [0.0066; 0.0116; 0.1116; 0.1296; 0.2166; 0.4376; 0.7006; 0.8146; 0.9766; 0.9946; 0.9956];

while ((max(abs(dR)) >= dRtol) && (iter < maxiter))
    z = zeros(length(displacement),1);
    for i = 1:length(z)
        z(i) = computeI(r,displacement(i));
    end
    
    e = VoltageFraction - z;
    
    J = jacobianS(@computeI,r,displacement);
    
    dR = J \ e;
    
    r = r + dR;
    
    iter = iter + 1;
    
end

if (iter > maxiter)
        disp('Reached max iterations!');
end

z = zeros(length(displacement),1);
for i = 1:length(z)
    z(i) = computeI(r,displacement(i));
end

e = VoltageFraction - z;
disp('magnitude of error = ');
error = norm(e)


estRbeam = r(1);
estRblock = r(2);
end

% computes the jacobian of the vector valued function at x0
function [J] = jacobianS(f, x0, D)

dx = 1e-7;

J = zeros(length(D),length(x0));

for j = 1:length(D)
    f0 = f(x0,D(j));
    for i = 1:length(x0)
        x = x0;
        % Perturb just the i'th value of x
        x(i) = x(i) + dx;
        
        J(j,i) = (f(x,D(j)) - f0) / dx; 
    end
end

end

function [I] = computeI(R, d)
    
    I = irradianceEstimate(R(1), R(2), d);
end

