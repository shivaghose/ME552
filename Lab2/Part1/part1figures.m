% generates the figures for part 1 that we care about

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compares the purely theoretical radiance fraction with the actual
% observed data

rbeam = 2.3650; % units in mm
rblock = 6.350;

% actual data follows
displacement = [3.64; 4.6; 5.07; 5.5; 5.55; 6; 6.5; 7; 7.45; 7.93; 9.36];
VoltageFraction = [0.0066; 0.0116; 0.1116; 0.1296; 0.2166; 0.4376; 0.7006; 0.8146; 0.9766; 0.9946; 0.9956];

irradianceSweep(rbeam, rblock);

hold on;
plot(displacement, VoltageFraction, 'k');
legend('Theoretical Ratio', 'Measured Ratio');
hold off;