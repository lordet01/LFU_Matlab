function bit_write(b_seq, fbit)

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

for k = 1:d_L
    if b_seq(k) == 1
        fprintf(fbit, '%c', '1');
    else
        fprintf(fbit, '%c', '0');
    end
end
