filenames = {'processed1PI.xls','processed1_5PI.xls','processed2PI.xls','processd3PI.xls','processed4PI.xls'};
startVal = [0, 0, 0, 0, 0];
endVal = [1*pi, 1.5*pi, 2*pi, 3*pi, 4*pi];
name = {'1 Pi','1.5 Pi','2 Pi','3 Pi','4 Pi'};

for f = 1:length(filenames)
    VAL = xlsread(filenames{f});
    T = VAL(:,2); 
    X = VAL(:,3);
    
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
    ylabel('Angular Position (radians)');
    title(strcat('Response of Position Controller for Step of ',name{f},' Radians'));
    hold off
end