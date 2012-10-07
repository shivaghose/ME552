function [z1plus, z2plus, z3plus] = stateEstimator(x, z1, z2, z3, h, a1, a2, b1, b2, b3, d)
% inspired by Y. X. Su's "A simple Nonlinear Velocity Estimator for
% High-Performance Motion Control
% x, z1, z2, z3, h, alpha1, alpha2, beta1, beta2, beta3, delta

% h is the sampling step
% z1, z2, z3 are the results of the previous computation
% delta, a1, a2, b1, b2, b3 are design parameters
% arbitrary parameters are arbitrary


e = z1 - x;
z1plus = z1 + h * (z2 - b1 * e);
z2plus = z2 + h * (z3 - b2 * fal(e, a1, d));
z3plus = z3 - h * b3 * fal(e, a2, d);

end

function [out] = fal(e, a, delta)
% inputs are epsilon, alpha, delta

if abs(e) > delta
    out = sign(e)*(abs(e))^a;
else
    out = e*(delta^(a - 1));
end
end