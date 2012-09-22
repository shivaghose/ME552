function [F] = estimateForceSimple(A,B,I,x)
% a much simpler model of force exterted by the electromagnet

F = -(A * I * I) / ((B + x) * (B + x));
end

