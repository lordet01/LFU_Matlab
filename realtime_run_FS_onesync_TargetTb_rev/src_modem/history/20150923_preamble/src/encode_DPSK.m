function [y, b_prev frame_idx]=encode_DPSK(b_seq, b_prev, frame_idx, fbit, w, Tb, Fs, fc, AMP)

[seq_L,d_L] = size(b_seq);
if seq_L > 1
    tmp = zeros(1,d_L*seq_L);
    tmp(1:d_L) = b_seq(1,:);
    for ii = 2: seq_L
        tmp(d_L*(ii-1)+1 : d_L*ii) = b_seq(ii,:);
    end
    b_seq = tmp;
    d_L = d_L * seq_L;
end

db = zeros(d_L,1);
y = zeros(1,d_L*Tb);
for k = 1:d_L
    if b_seq(k) == 1
        fprintf(fbit, '%c', '1');
    else
        fprintf(fbit, '%c', '0');
    end
    
    db(k) = ~xor(b_seq(k), b_prev);
    
%     t_shift = frame_idx * Tb;
    t_shift = 0;
    if db(k) == 1
        s = (AMP * w.*cos(2*pi*fc*(1+t_shift:Tb+t_shift)/(Fs)));
    else
        s = -(AMP * w.*cos(2*pi*fc*(1+t_shift:Tb+t_shift)/(Fs)));
    end
    

    y((k-1)*Tb+1 : k*Tb) = s;

    frame_idx = frame_idx+1;
    b_prev = db(k);
end