% generates the figures for part 1 that we care about

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compares the purely theoretical radiance fraction with the actual
% observed data

rbeam = 2.3650; % units in mm
rbeam = 1.6706; % only fitting estimateRbeam
rblock = 6.350;

rbeam = 1.398935; %fitting both
rblock = 6.070821; 

% actual data follows
displacement = [3.49; 4.45; 4.927; 5.35; 5.4; 5.85; 6.35; 6.85; 7.3; 7.78; 9.21];
VoltageFraction = [0.0066; 0.0116; 0.1116; 0.1296; 0.2166; 0.4376; 0.7006; 0.8146; 0.9766; 0.9946; 0.9956];

irradianceSweep(rbeam, rblock);

hold on;
plot(displacement, VoltageFraction, 'k');
legend('Theoretical Ratio', 'Measured Ratio');
hold off;