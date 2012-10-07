function [] = optimizeEstimator()
% for optimizing the parameters to the state estimator

% first going to make sure that things actually work
samplingPeriod = 1 / 200; % sampling at 200 Hz
resolution = 2000; % ticks per revolution

T = [0 : samplingPeriod : 5]; % time

[Ft, dFt] = Sine(T); % actual values
encoderFt = apply(@discretize, @Sine, T, resolution); % encoder reading
% figure
% hold on;
% plot(T,Ft,'k')
% plot(T,encoderFt,'b')
% legend('Actual Position','Reported Encoder Position')
% xlabel('Time (seconds)')
% ylabel('Position (radians)')
% title('Encoder Emulation')
% hold off;

%%%%%%%%%%%%%%%%%%%%%
% Optimization
%%%%%%%%%%%%%%%%%%%%%

% our initial guess for parameters
% a1, a2, b1, b2, b3, d
h = 0.0001;
x = [0.5; 0.25; 0.5; 0.5; 0.5; 0.5];

dxtol = 1e-6;
dX = inf;
iter = 0;
maxiter = 1000;

while ((max(abs(dX)) >= dxtol) && (iter < maxiter))
    
    % compute errors
    [Z1, Z2] = shootEstimator(T, encoderFt, h, x);
    E = computeObsError(Ft, dFt, Z1, Z2);
    
    disp('mag E = ')
    norm(E)
    size(E)
    
    % compute Jacobian matrix
    J = computeJacobian(T, encoderFt, h, x);
    size(J)
    
    dX = J \ E
    
    x = x + dX;
    
    iter = iter + 1;
end

[Z1, Z2] = shootEstimator(T, encoderFt, h, x);

disp('magnitude of final error = ')
norm(computeObsError(Ft,dFt, Z1, Z2));


%%%%%%%%%%%%%%%%%%%%%%%
% Make some plots
%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on;
plot(T,Ft,'k')
plot(T,Z1,'b')
legend('Actual Position', 'Filtered Encoder Positon')
xlabel('Time (seconds)')
ylabel('Position (radians)')
hold off

figure
hold on;
plot(T, dFt,'k');
plot(T, Z2, 'b');
legend('Actual Velocity', 'Filtered Encoder Velocity')
xlabel('Time (seconds)')
ylabel('Angular Velocity (rad/s)')
hold off

end

function [J] = computeJacobian(T, Ft, h, x)

J = zeros(2 * length(T), length(x));

delta = 1e-6;

[Z10 Z20] = shootEstimator(T, Ft, h, x);

% compute nominal result
pZ1 = zeros(length(T),length(x));
pZ2 = zeros(length(T),length(x));

% compute perturbed result
for v = 1 : length(x)
    px = x;
    px(v) = px(v) + delta;
    [pZ1(:,v) pZ2(:,v)] = shootEstimator(T, Ft, h, px);
end

% compute partial derivative 
for i = 1 : length(T)
    for v = 1 : length(x)
        J(2 * i - 1, v) = (Z10(i) - pZ1(i,v)) / delta;
        J(2 * i, v) = (Z20(i) - pZ2(i,v)) / delta;
    end
end

end

function [E] = computeObsError(Ft, dFt, Z1, Z2)
% have to interleave the values for the error computation
assert((length(Ft) == length(dFt)) && (length(Ft) == length(Z1)) && (length(Ft) == length(Z2)), 'Should have been the same length'); 
E = zeros(2*length(Ft),1);
for i = 1 : length(Ft)
    E(2 * i - 1) = Ft(i) - Z1(i);
    E(2 * i) = dFt(i) - Z2(i);
end
end

function [Z1 Z2] = shootEstimator(T, F, h, x)
% extract variables
a1 = x(1); a2 = x(2); b1 = x(3); b2 = x(4); b3 = x(5); d = x(6);

Z1 = zeros(length(T),1);
Z2 = zeros(length(T),1);

z1minus = 0;
z2minus = 0;
z3minus = 0;
for i = 1 : length(T)
    [z1plus, z2plus, z3plus] = stateEstimator(F(i),z1minus,z2minus,z3minus,h,a1,a2,b1,b2,b3,d);
    Z1(i) = z1plus;
    Z2(i) = z2plus;
    % save away these values for the next round
    z1minus = z1plus; z2minus = z2plus; z3minus = z3plus;
end

end

function [out] = apply(fun, target, T, resolution)
out = zeros(length(T),1);
for i = 1 : length(T)
    out(i) = fun(target, T(i), resolution);
end

end

function [x xdot] = Sine(t)
x = sin(t);
xdot = cos(t);
end