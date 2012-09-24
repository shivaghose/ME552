% parameters for magnet
%A = 8.9574e-06;
%B = 1.8741e-05;
Ac = 1.9752e-05;
Bc = 8.9294e-04;
A = 2.2795e-05;
B = 0.0011656;

% good values
%A = 3.4963e-04; B = 0.0011;

% good data set
%I = [0.358, 0.308, 0.365, 0.298, 0.288, 0.277] ;
%X = [3.23, 2.73, 3.235, 3.13, 3.035, 2.935] ./ 1000;

Ishiva = [0.308, 0.277, 0.288, 0.298, 0.358, 0.365];
Xshiva = [2.73, 2.935, 3.035, 3.135, 3.23, 3.235] ./ 1000;

X = [0.52; 1.95; 2.38; 2.88; 3.38; 3.88; 4.31] ./ 1000; % want m
V = [0.74; 1.17; 1.349; 1.565; 1.74; 2; 2.25];

% compute current
I = actualdriverVtoI(V);

x = (0:0.01:4.5);
x = x ./ 1000;
f = -0.0981;
Iest = zeros(length(x),1);
for i = 1:length(x)
    Iest(i) = requiredCurrent(A,B,x(i));
end

ICest = zeros(length(x),1);
for i = 1:length(x)
    ICest(i) = requiredCurrent(Ac,Bc,x(i));
end

figure;
hold on;
plot(X,I,'ko');
plot(x,Iest,'b');
plot(x,ICest,'r');
%plot(Xshiva,Ishiva,'r');
title('Necessary Magnet Current: Experimental vs Fit Values');
xlabel('X (m)');
ylabel('Current (A)');
legend('Observed Position & Current','Necessary Current: Force Fit', 'Neccessary Current: Current Fit')
hold off;



%%%%%%%%%%%
% Driver Circuit Data

Vin = [0; 0.5; 1; 1.25; 1.5; 1.75; 1.8; 1.85; 1.9; 1.95; 2; 2.05; 2.1; 2.15; 2.2; 2.25; 2.3; 2.4; 2.45; 2.5; 2.6; 2.7];
Iout = [-0.015; 0.071; 0.159; 0.203; 0.248; 0.292; 0.302; 0.311; 0.32; 0.329; 0.338; 0.348; 0.357; 0.366; 0.375; 0.384; 0.394; 0.401; 0.4; 0.404; 0.4; 0.399];

Vinm = 0:0.1:3;
Ioutm = 0.1781.*Vinm - 0.0179;
Ioutt = 0.042331.*Vinm;

figure
hold on;
plot(Vin,Iout,'ko')
plot(Vinm,Ioutm,'b')
plot(Vinm,Ioutt,'r')
legend('Measured Ouput Currents','Fit Current Model','Theoretical Current Model')
title('Driver Circuit Modeling')
xlabel('Input Voltage (V)')
ylabel('Electromagnet Current (I)');
hold off;
