function [Pest] = fitexperimental1()
% uses parameterEstimation matlab function to fit the data

%P0 = [8.5e-6; 3.7e-6; (5.6e-3)/2; 4.24e-2]; % initial guesses
P0 = [0.0001; 0.0001; 0.1192; 0.0467] % better initial guess

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

Pest = parameterSweep(T',Theta',I,P0)
%Pest = parameterEstimation(0,T',Theta',I,P0)

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