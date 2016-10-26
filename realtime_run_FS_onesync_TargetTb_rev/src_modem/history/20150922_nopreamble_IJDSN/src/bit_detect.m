function [bit] = bit_detect(data_prev, data_now)

THR = 0;
th = zeros(1, 8);

h = data_prev .* data_now;

th = sum(h);
if th > THR
    bit = 1;
else
    bit = 0;
end

