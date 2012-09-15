function [] = irradianceSweep(rbeam, rblock)
% sweeps the distances as a function of beam and ball radius to generate a 
% pretty plot

assert((rbeam >= 0) && (rblock >= 0), 'The radii should be non-negative');

num = 1000;

dd = (rbeam + rblock) / (num - 1);

D = 0: dd: (rbeam + rblock);
IR = zeros(1,num);

length(D)
length(IR)

for i = 1: num;
    IR(i) = irradianceEstimate(rbeam,rblock,D(i));
end

figure;
plot(D,IR);
xlabel('Distance between centers');
ylabel('Fraction of Maximum Irradiance');
title('Irradiance Curve');

end

