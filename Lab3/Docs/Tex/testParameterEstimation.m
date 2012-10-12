function [Pest] = testParameterEstimation()
% uses parameterEstimation matlab function to fit the data

J = 0.000041; % initial guesses
Kt = 0.042;
Bforwards = 5.4288e-4 * Kt; % 10*5.4288e-4 * Kt
Breverse = 6.728054e-4 * Kt; % 0.9*6.728054e-4 * Kt
Tfforwards = 1.75*0.4637 * Kt; %1.75*0.4637 * Kt
Tfreverse = 0.6071722 * Kt; % 

% Bforwards = 0.000023725941605;
% Breverse = 0.000020265029631;
% Tfforwards = 0.034057546192782;
% Tfreverse = 0.025282905174523;
% 
% J = 0.000043471112655;
% Kt =  0.042000000000000;
% Bforwards =  0.000002; % negative values are no good -0.000011421939988
% Breverse = 0.000002;
% Tfforwards = 0.033926022489524;
% Tfreverse = 0.025288849161596;
% 
% J = 0.000046420828236;
% Kt =    0.044062508073865;
% Bforwards =    0.000001678546329;
% Breverse =   0.000002112293049;
% Tfforwards =   0.034540796108616;
% Tfreverse =   0.026025786119668;

% best so far with magnitude of error of 40.2317373
%  0.000044142993118
%   0.042016985993061
%  -0.000018376326419
%   0.000000158849427
%   0.033878935596711
%   0.025319881527519

% best numbers!!! magnitude of error of 27.649 on short data set with
% 2*10^4 points
J = 0.000035780751317;
Kt = 0.032448801388174;
Bforwards = 0.000032339885634;
Breverse = 0.000041519670783;
Tfforwards = 0.022294359068982;
Tfreverse = 0.015632909278982;

% error of 23.4941109 on short data set with 2*10^4 points
% 0.000035051428588
%    0.031449871029185
%    0.000010858577896
%    0.000048592997081
%    0.022438919281204
%    0.014309603079420

% initial error with best numbers above 8.425662612875667e+02
% on medium came up with this 3.259151910806145e+02
% 0.000028710305642
%    0.031941181768441
%    0.000075360245042
%    0.000024930146369
%    0.019532312877692
%    0.016530560525051

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

J = 6.06724e-5;
B = 4.05562e-5;
Kt = 0.01902;
Bforwards = B;
Breverse = B;
Tf = 0.000145312;
Tfforwards = Tf;
Tfreverse = Tf;

%Pest = parameterSweep(T',Theta',I,P0)
%Pest = nonuniformParameterEstimation(0,T,Theta,I,P0,);
state0 = [J; Kt; Bforwards; Breverse; Tfforwards; Tfreverse];
locked = [ 1; 1;     1    ;     1   ;     1    ;       1     ];
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