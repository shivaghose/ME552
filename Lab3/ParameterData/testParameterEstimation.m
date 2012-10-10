function [Pest] = testParameterEstimation()
% uses parameterEstimation matlab function to fit the data

J = 0.000041; % initial guesses
Kt = 0.042;
Bforwards = 5.4288e-4 * Kt; % 10*5.4288e-4 * Kt
Breverse = 6.728054e-4 * Kt; % 0.9*6.728054e-4 * Kt
Tfforwards = 1.75*0.4637 * Kt; %1.75*0.4637 * Kt
Tfreverse = 0.6071722 * Kt; % 



VAL = xlsread('shortSinusoidalInput.xls');

T = VAL(:,2); % delta Time
Theta = VAL(:,3); % need negative
Voltage = VAL(:,4);

% current is related to voltage
I = Voltage;

% % plot the experimental data
% figure 
% hold on
% plot(T,Theta,'b');
% plot(T,I,'r');
% xlabel('Time (seconds)');
% legend('Angular Position','Motor Current');
% hold off

%Pest = parameterSweep(T',Theta',I,P0)
%Pest = nonuniformParameterEstimation(0,T,Theta,I,P0,);
state0 = [J; Kt; Bforwards; Breverse; Tfforwards; Tfreverse];
locked = [ 1; 1;     0    ;     0   ;     0    ;       0     ];
Pest = enchancedParameterEstimation(T,Theta,I,state0,locked)
end

function [X, L] = trimData(X)
L = length(X);
for i = 1:length(X)
    if isnan(X(i))
        L = i - 1;
        break;
    end
end
X = X(1:L);
end