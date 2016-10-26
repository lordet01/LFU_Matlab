function bit_blk_c = matB2C(bit_blk, M, C)

[N, B] = size(bit_blk);

seq = zeros(1, N*B);
for ii = 1:N
   seq((ii-1)*B+1 : ii*B) = bit_blk(ii,:); 
end

bit_blk_c = zeros(M,C);
for jj = 1:M
    bit_blk_c(jj,:) = seq((jj-1)*C+1 : jj*C);
end
