function [fullR] = computeR(Rf_, Rg_, a)
% for Lab part 2, given Rf, and Rg, and a = [a1, a2, a3, a4]
% note that all a1 ... a4 must be positive
% compute R1, R2, R3, R4

% just using numerical derivative

    global Rf Rg;
    Rf = Rf_; Rg = Rg_;
    
    % do nonlinear least squares
    dRtol = 1e-9;
    maxiter = 1000;

    % construct R
    R = [1; 1; 1; 1]; % 1 as initial estimate for each
    
    dR = Inf; % initialize so that we enter the loop
    iter = 0;
    
    while (max(abs(dR)) >= dRtol && iter <= maxiter)
        
        % compute error
        e = a - computeA(R);
        
        %error = norm(e)
        
        % compute Jacobian at the current estimate
        J = jacobian(@computeA, R);
    
        % compute update matlab handles the system solving
        dR = J \ e;
    
        % update
        R = R + dR;
    
        iter = iter + 1;
    
    end
    
    if (iter > maxiter)
        disp('Reached max iterations!');
    end
    
    disp('A = ')
    computeA(R)
    
    e = a - computeA(R);
    disp('magnitude of error');
    error = norm(e)
    
    % now assemble the full R
    fullR = [R; Rf; Rg];
end

% where R = [R1, R2, R3, R4]

% computes the jacobian of the vector valued function at x0
function [J] = jacobian(f, x0)

    dx = 1e-7;
    
    f0 = f(x0);
    
    for i = 1:length(x0)
        x = x0;
        % Perturb just the i'th value of x
        x(i) = x(i) + dx;     
        
        df(:,i) = f(x) - f0;  % delta f
    end

    J = df/dx;
end

function [a] = computeA(R) 
    a = [computeA1(R); computeA2(R); computeA3(R); computeA4(R)];
end

function [a1] = computeA1(R) 
    % extract values
    R1 = R(1); R2 = R(2); R3 = R(3); R4 = R(4); 
    global Rf Rg;
    
    % compute a1
    a1 = (Rf / R3 + Rf / R4 + 1) * 1 / (R1 / Rg + R1 / R2 + 1);
end

function [a2] = computeA2(R) 
    % extract values
    R1 = R(1); R2 = R(2); R3 = R(3); R4 = R(4); 
    global Rf Rg;
    
    a2 = (Rf / R3 + Rf / R4 + 1) * 1 / (R2 / R1 + R2 / Rg + 1);
    
end

function [a3] = computeA3(R)
    % extract values
    R3 = R(3); 
    global Rf;
    
    a3 = Rf / R3;
end

function [a4] = computeA4(R)
    % extract values
    R4 = R(4); 
    global Rf;
    
    a4 = Rf / R4;
end


