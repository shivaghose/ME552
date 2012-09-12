% generates the expected and actual bode plots of our system for a variety
% of frequencies

% notes on compensator.vi
% period = 1ms
% input amplitude ~ 0.1 depending on frequency
% multiplying system output by -0.01 before storing it because system has a
% very large amplification

% phase extracted using 2 tone blocks to measure the phase and then taking
% the difference of these extracted phases.  becomes unstable at high
% frequencies because can no longer fully reconstruct the wave due to
% limited sampling frequency

% using tone measurements to extract amplitude of waves and then dividing
% these to find the ratio of input to output

% experimental data
% frequecy in Hz
F = [0.01; 0.05; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0; 
    2.0; 5.0; 8.0; 10; 15; 20; 30; 40; 50; 60; 70; 80; 90; 100; 200; 300];
F = F*2*pi;

exPhase = [-19.301; -5.610; -1.375; 7.031; 8.116; 11.514; 14.619; 18.203; 
    20.195; 22.838; 25.151; 27.405; 42.054; 50.611; 46.378; 41.544; 
    31.735; 21.8406; 8.479; 2.399; -4.91; -11.186; -16.574; 
    -21.736; -26.563; -30.648; -30; -30]; % extra phase data is messing us up286.687; -101.823

AmplitudeRatio = [0.817; 0.583; 0.572; 0.571; 0.577; 0.585; 
    0.597; 0.608; 0.622; 0.637; 0.654; 0.675; 0.908; 1.652; 2.332; 
    2.451; 3.067; 3.144; 3.194; 4.062; 4.083; 4.099; 4.111; 4.112; 
    4.120; 4.140; 5.07625; 5.04914];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the real transfer function
s = tf('s');
G = ((1 + 0.1*s) * (1 + 5*s)) / ((1 + 0.01*s) * (1 + 10*s));

[Mag, Phase, W] = bode(G, {0.001, 500*2*pi});
MdB = zeros(length(Mag),1);
P = zeros(length(Phase),1);
for i = 1:length(Mag)
    MdB(i,1) = 20*log10(Mag(1,1,i));
    P(i,1) = Phase(1,1,i);
end

% convert to dB
aRdB = 20*log10(AmplitudeRatio);

% h = gcf; % get handle to current figure
% axesObjs = get(h, 'Children')
% dataObjs = get(axesObjs, 'Children')
% 
% % 2 groups of data one for each figure
% magLine = get(dataObjs{1}(1),'Children');
% magX = get(magLine, 'XData');
% magY = get(magLine, 'YData');
% 
% phaseLine = get(dataObjs{3}(1),'Children');
% phaseX = get(phaseLine,'XData');
% phaseY = get(phaseLine,'YData');

% now plot our data
figure;
handle = subplot(2,1,1);
get(handle)
hold on
plot(F,aRdB);
plot(W,MdB,'--');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
set(handle,'XScale','log');
hold off;

handle = subplot(2,1,2);
hold on
semilogx(F,exPhase);
semilogx(W,P,'--');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
set(handle,'XScale','log');
hold off;
