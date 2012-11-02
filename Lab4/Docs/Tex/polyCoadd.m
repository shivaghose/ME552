function [pSum] = polyCoadd(p1, p2)
% returns sum of two arrays representing polynomial coefficients

if length(p1) >= length(p2)
    pSum = p1;
    for i = length(p2):-1:1
        pSum(length(pSum) - i + 1) = pSum(length(pSum) - i + 1) + p2(length(p2) - i + 1);
    end
else
    pSum = p2;
    for i = length(p1):-1:1
        pSum(length(pSum) - i + 1) = pSum(length(pSum) - i + 1) + p1(length(p1) - i + 1);
    end
    
end
end