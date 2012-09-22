% parameters for magnet
A = 8.9574e-06;
B = 1.8741e-05;

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
I = driverVtoI(V);

x = (0:0.01:4);
x = x ./ 1000;
f = -0.0981;
Iest = zeros(length(x),1);
for i = 1:length(x)
    Iest(i) = requiredCurrent(A,B,x(i));
end

figure;
hold on;
plot(X,I,'k');
plot(x,Iest,'b');
plot(Xshiva,Ishiva,'r');
title('Necessary Magnet Current: Experimental vs Fit Values');
xlabel('X (m)');
ylabel('Current (A)');
hold off;
