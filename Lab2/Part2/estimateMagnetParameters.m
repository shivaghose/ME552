function [A,B] = estimateMagnetParameters()
% using data taken from the magnet, tries to fit the simple force
% model to the data

% actual data
X = [0.52; 1.95; 2.38; 2.88; 3.38; 3.88; 4.31] ./ 1000; % want m
V = [0.74; 1.17; 1.349; 1.565; 1.74; 2; 2.25];

% compute current
I = actualdriverVtoI(V);


% different data set
%I = [0.358, 0.308, 0.365, 0.298, 0.288, 0.277] ;
%X = [3.23, 2.73, 3.235, 3.13, 3.035, 2.935] ./ 1000;

% Aplus = 8.9387e-9 Bplus = -0.0031
% Aminus = 7.6128e-7 Bminus = -0.0022
% note that using Aminus and Bminus above yielded a lower residual
% 7.6128e-7; 1.8741e-5
x = [1.97e-6; 8.929e-4]; % inital estimate of A and B;
% actual A = 8.9574e-06 B = 1.8741e-5

dxtol = 1e-11;
maxiter = 1000;

dx = Inf;
iter = 0;

% actual observations
f =  - 0.01 * 9.81; % actual F is the same for each
% negative to be consistent with sign system
Fobs = ones(length(X),1).*f;

while ((max(abs(dx)) >= dxtol) && (iter < maxiter))
    % compute force errors
    F = zeros(length(X),1);
    for i = 1: length(F)
        F(i) = estimateForceSimple(x(1), x(2), I(i), X(i)); % A, B, I, distance
    end
    E = Fobs - F;
    
    magE = norm(E);
    
    % compute jacobian
    J = jacobianS(x, I, X);
    
    dx = J \ E;
    
    x = x + dx;
    magD = norm(dx)
    
    iter = iter + 1;
    
end

if (iter > maxiter)
    disp('Reached max iterations!');
end

% compute force errors
F = zeros(length(X),1);
for i = 1: length(F)
    F(i) = estimateForceSimple(x(1), x(2), I(i), X(i)); % A, B, I, distance
end
E = Fobs - F;

disp('magnitude of final error = ');
error = norm(E)

A = x(1);
B = x(2);

end

% computes a special jacobian, only optimizing A and B, assuming
% observations are perfect
function [J] = jacobianS(x0, I, X)

dx = 1e-8;

% partial with respect to A and then B
J = zeros(length(I),length(x0));
for i = 1:length(x0) % column
    
    % perterbation
    deltax = zeros(length(x0),1);
    deltax(i) = dx;
    
    % perterbed x
    xp = x0 + deltax;
    
    for j = 1: length(I) % row
        % force at current A and B
        f0 = estimateForceSimple(x0(1),x0(2),I(j),X(j));
        
        % perterbed force
        fp = estimateForceSimple(xp(1),xp(2),I(j),X(j));
        
        J(j,i) = (fp - f0) / dx;
        
    end
end
J
end

