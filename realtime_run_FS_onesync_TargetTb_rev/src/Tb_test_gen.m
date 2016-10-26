function Tb_test_gen(snd_path, rcv_path, snd_list, A_prm, Tb, p)

%%Parameters
fs = p.fs;
fc = p.fc;
g = p.g; %Spreading factor
% A0 = p.A0; %Unit-normalizing constant for Attenuation factor
l = 0.0;

str =[];

%Load Tx files
for i = 1:length(snd_list)
    fname_snd = snd_list(i).name;
    s_full = wavread([snd_path,'/',fname_snd]);
    s_full = s_full .* 32767;
    s = [zeros(p.sil_len,1); s_full; zeros(p.sil_len,1)];
%     s_full = s_full(p.comp_len+1:2*p.comp_len);
	K = length(s);
    K21 = K/2+1;
    
    %% Attenuation function
    fc_bin = fc / fs * 2 * K21;
    funA = @(x,xdata)(x(1) * xdata.^g) .* exp(xdata .* (x(2) * fc_bin .^ x(3)));
    A = funA(A_prm, l);

    %Notify chosen Tb
    Tb_out = Tb(i);
    
    % Generate modeld Rx signals
    l_Pa = 10*log10(sum(A)+1);
    Pa = 10^(l_Pa / 10);
    s_a = s ./ sqrt(Pa); %Attenuated Tx signal by distance
    rcv_sig = s_a;
    
    %Add dummy silences
    rcv_sig_out = [zeros(p.N_guard_len,1); rcv_sig; zeros(p.N_guard_len,1)];
    rcv_sig_out = rcv_sig_out  ./ 32767;
    
    fname_rcv = snd_list(i).name;
    fname_rcv = strtok(fname_rcv,'.');
    fname_rcv_out = [rcv_path,'/Tb_Test/',fname_rcv,'.wav'];
    disp(fname_rcv_out);
    wavwrite(rcv_sig_out, p.fs, fname_rcv_out);
    
    fname_txt_rcv = ['src_modem/output/ref/',fname_rcv,'.txt'];
    bit_name_txt = ['src_modem/bit/ref/',fname_rcv,'.txt'];
    preamble_rcv = ['src_modem/preamble/p_',fname_rcv,'.pcm.mat'];
    u_demod(fname_rcv_out, fname_txt_rcv, bit_name_txt, preamble_rcv, Tb_out, p);
end

end
