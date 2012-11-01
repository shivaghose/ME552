function [Pest] = fitexperimental1()
% uses parameterEstimation matlab function to fit the data

P0 = [8.5e-6; 3.7e-6; 0; 4.24e-2]; % initial guesses
%P0 = [0.0001; 0.0001; 0.1192; 0.0467]; % better initial guess
%P0 = [0.00015; 0.00010; 0.001; 0.0467];
%P0 = [0.000303362119919; 0.000202781054320; 0.0243853568; 0.09510344];
%P0 = [0.000060672423984; 0.000040556210864; 0.014631214080000; 0.019020688000000];
%P0 = [1e-5; 1e-5; 1e-4; 0.04];

VAL = xlsread('ExperimentalData1.xls');

T = VAL(:,2);
Theta = VAL(:,3);
AngularVelocity = VAL(:,4);
Voltage = VAL(:,5);

% trim data
[AngularVelocity L] = trimData(AngularVelocity);
[Voltage L] = trimData(Voltage);
T = T(1:L)';
Theta = Theta(1:L)';

% current is related to voltage
I = Voltage;


% plot the experimental data
figure 
hold on
plot(T,Theta,'b');
plot(T,I,'r');
xlabel('Time (seconds)');
legend('Angular Position','Motor Current');
hold off

size(T')
size(Theta')
size(I)

%Pest = parameterSweep(T',Theta',I,P0)
Pest = parameterEstimation(0,T',Theta',I,P0)

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