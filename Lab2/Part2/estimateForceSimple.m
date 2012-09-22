function [F] = estimateForceSimple(A,B,I,x)
% a much simpler model of force exterted by the electromagnet
% correct values for A = 8.9574e-6  B = 1.8741e-5
% from regression on the floating magnet current measurements ||E|| =
% 0.0396

F = -(A * I * I) / ((B + x) * (B + x));
end

