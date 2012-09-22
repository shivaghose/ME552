function [A, B] = estimateMagnetParametersBetter()
% instead of trying to minimize error in the force domain, minimize in
% current

% these were lift voltages vs height
X = [0.52; 1.95; 2.38; 2.88; 3.38; 3.88; 4.31] ./ 1000; % want m
V = [0.74; 1.17; 1.349; 1.565; 1.74; 2; 2.25];

X
% compute current
I = driverVtoI(V);

global f;
f = - 0.01 * 9.81; % the weight of the sphere

x = [1e-6; 0.01]; % initial guess for A and B

dxtol = 1e-10;
maxiter = 1000;

dx = Inf;
iter = 0;

Iest = zeros(length(I),1);
while ((max(abs(dx)) >= dxtol) && (iter < maxiter))
    % compute current error
    for i = 1:length(I)
        Iest(i) = requiredCurrent(x(1),x(2),X(i));
    end
    E = I - Iest;
    
    J = jacobianS(x,X);
    
    dx = J \ E;
    
    x = x + dx;
    magD = norm(dx)
    
    iter = iter + 1;
end

if (iter > maxiter)
    disp('Reached max iterations!');
end

% last error check
% compute current error
for i = 1:length(I)
    Iest(i) = requiredCurrent(x(1),x(2),X(i));
end
E = I - Iest;
disp('magnitude of final error = ')
magE = norm(E)

A = x(1)
B = x(2)

end

function [J] = jacobianS(x0, X)
dx = 1e-8;
J = zeros(length(X),length(x0));

for i = 1:length(x0) % column
    % perterbation
    % perterbation
    xp = x0;
    xp(i) = xp(i) + dx;
    
    for j = 1: length(X) % row
        % current
        i0 = requiredCurrent(x0(1),x0(2),X(j));
        % perturbed current
        ip = requiredCurrent(xp(1),xp(2),X(j));
        
        J(j,i) = (ip - i0) / dx;
    end
end

end

function [I] = requiredCurrent(A,B,x)
% how much current do we need to hold it up at this distance according to
% this model
global f;
I = sqrt(-f * (B+x)*(B+x) / A);

end

