function [estRbeam] = estimateRbeam()
% a bit of non-linear least squares to choose rbeam that best explains the
% data that we've seen

global rBlock;

rBlock = 6.35;

% initial estimate for R = [rbeam];
r = [2.3650];

dRtol = 1e-9;
maxiter = 1000;

dR = Inf;
iter = 0;

% our observed data
displacement = [3.64; 4.6; 5.07; 5.5; 5.55; 6; 6.5; 7; 7.45; 7.93; 9.36];
VoltageFraction = [0.0066; 0.0116; 0.1116; 0.1296; 0.2166; 0.4376; 0.7006; 0.8146; 0.9766; 0.9946; 0.9956];

while (max(abs(dR) >= dRtol) && (iter < maxiter))
    z = zeros(length(displacement),1);
    for i = 1:length(z)
        z(i) = computeI(r(1),displacement(i));
    end
    
    e = VoltageFraction - z;
    
    J = zeros(length(displacement),1);
    for i = 1:length(J)
        J(i) = jacobianS(@computeI, r(1), displacement(i));
    end
    
    dR = J \ e;
    
    r = r + dR;
    
    iter = iter + 1;
    
end

if (iter > maxiter)
        disp('Reached max iterations!');
end

z = zeros(length(displacement),1);
for i = 1:length(z)
    z(i) = computeI(r(1),displacement(i));
end

e = VoltageFraction - z;
disp('magnitude of error = ');
error = norm(e)


estRbeam = r(1);
end

% computes the jacobian of the vector valued function at x0
function [J] = jacobianS(f, x0, d)

    dx = 1e-7;
    
    f0 = f(x0,d);
    
    for i = 1:length(x0)
        x = x0;
        % Perturb just the i'th value of x
        x(i) = x(i) + dx;     
        
        df(:,i) = f(x,d) - f0;  % delta f
    end

    J = df/dx;
end

function [I] = computeI(R, d)
    global rBlock;
    
    I = irradianceEstimate(R(1), rBlock, d);
end