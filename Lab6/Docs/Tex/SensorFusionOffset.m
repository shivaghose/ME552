function [ offset, start] = SensorFusionOffset( Encoder, filteredData )
% given encoder measurements, figures out where the first sign of motion is
% and determines the zero rate in Volts
% gives ZeroRate in V and variance in V^2
% in this case we are looking for exactly 0

start = 1;

while (start < length(Encoder)) && (abs(Encoder(start)) == 0)
    start = start + 1;
end

offset = mean(filteredData(1:start));
end

