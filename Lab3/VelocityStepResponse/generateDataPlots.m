filenames = {'step0_5Pi.xls','stepPi.xls','step1_5Pi.xls','step2Pi.xls','step5Pi.xls','step10Pi.xls','step2PiNoise.xls'};
startVal = [0, 0, 0, 0, 0, 0, 0];
endVal = [0.5*pi, 1*pi, 1.5*pi, 2*pi, 5*pi, 10*pi, 2*pi];
name = {' 0.5 Pi', ' 1 Pi',' 1.5 Pi',' 2 Pi',' 5 Pi',' 10 Pi', ' 2 Pi with Noise'};

for f = 1:length(filenames)
    VAL = xlsread(filenames{f});
    T = VAL(:,2); 
    X = VAL(:,1);
    
    %%%%%%%%%%%%%%%%
    % Extract Details 
    [riseTime, settlingTime, percentOvershoot, SteadyStateErrorPercent] = systemDiagnosisStep(startVal(f),endVal(f),T,X);
    
    name(f)
    riseTime
    settlingTime
    percentOvershoot
    SteadyStateErrorPercent
    
    %%%%%%%%%%%%%%%%%
    % generate plot of step and response
    Step = endVal(f) .* ones(1,length(T));
    Step(1) = 0;
    
    figure
    hold on
    plot(T,Step,'k');
    plot(T,X,'b');
    xlabel('Time (seconds)');
    ylabel('Angular Velocity (radians/second)');
    title(strcat('Response of Velocity Controller for Step of ',name{f},' Radians'));
    hold off
end