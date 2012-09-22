function [Aplus Bplus Aminus Bminus] = estimateMagnetParametersAlt(I1,x1,I2,x2)

R = (I2 * I2) / (I1 * I1);

a = 1 - R;
b = 2 * (x2 - R * x1);
c = x2*x2 - R*x1*x1;

Bplus = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
Bminus = (-b - sqrt(b*b - 4*a*c)) / (2 * a);

f = 0.01 * 9.81;

Aplus = (f * (Bplus + x1)^2) / (I1 * I1);
Aminus = (f * (Bminus + x1)^2) / (I1 * I1);

end

