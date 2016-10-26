function [EbN0, B, C, SNR_cand, fname_rcv]=channel_sim_fs(snd_path, rcv_path, fname_snd, n, l, Tb_target, A_prm, fc, p)

%%Parameters
fs = p.fs;
g = p.g; %Spreading factor
% As = p.As; %Unit-normalizing constant for Attenuation factor


str =[];

%Load Tx files
% for i = 1:length(snd_name)
%     fname_snd = snd_list(i).name;
s_full = wavread([snd_path,'/',fname_snd]);
s_full = s_full .* 32767;
if length(s_full) >= p.analysis_len
    s = s_full(1:p.analysis_len);
else
    s = s_full;
end
s = [zeros(p.sil_len,1); s; zeros(p.sil_len,1)];
sample_num = length(s);

fc_L = fc - fs/Tb_target;
fc_H = fc + fs/Tb_target;


%% Get Eb/N0 for given l and SNRs
%Get PSD of Tx and Noise
S1 = abs(fft(s));
S = S1(1:sample_num/2+1).^2; %./ (sample_num*fs);
S(2:end-1) = 2*S(2:end-1);

N1 = abs(fft(n(1:sample_num)));
N = N1(1:sample_num/2+1).^2; %./ (sample_num*fs);
N(2:end-1) = 2*N(2:end-1);

%Get log N and S of subband
% l_S = 10*log10(S);
% l_N = 10*log10(N);


%Calculate band of interest
B_init_bin =  (1:sample_num*0.5+1);
% B_init_f = B_init_bin / (sample_num/2+1) * fs/2;

%Attenuation funfction
% % alpha = B_init_f .^ g;
% % A = As * l^g * exp(l * alpha);
% % A = A';
K = length(s_full);
K21 = K/2+1;

%% Attenuation function
fc_bin = fc / fs * 2 * K21;
funA = @(x,xdata)(x(1) * xdata.^g) .* exp(xdata .* (x(2) * fc_bin .^ x(3)));
A = funA(A_prm, l);

    
    

f_bin_L = round(2/fs * (sample_num/2+1) * fc_L);
f_bin_H = round(2/fs * (sample_num/2+1) * fc_H);

if f_bin_H > length(B_init_bin);
    BW_R = length(B_init_bin);
else
    BW_R = f_bin_H;
end
BW_L = f_bin_L;

B_cand =  B_init_bin(BW_L : BW_R);
SNR_cand = 10*log10(sum(S(B_cand))) - 10*log10(sum(N(B_cand))) - 10*log10(sum(A)+1);
B_width = B_cand;

% Calculate Capacity Eq.(18)
AN = (A+1) .* N(B_width);
K = S(B_width) + AN;
C = sum(log2( K./( AN )));
penalty_mod = p.BPS_num ./  Tb_target;
C = C .* penalty_mod ./ p.C_scale;

% Calculate Eb/N0 Eq.(24)
B = length(B_width) / (sample_num/2+1) * fs/2;
% CB = C / B;
EbN0 = 10*log10(B) -  10*log10(C) + SNR_cand;

% %Notify chosen Tb
% Tb_idx = i;

%Display
B = length(B_width) / (sample_num/2+1) * fs/2;
fprintf(repmat('\b',1,length(str)));
str = sprintf('Tb %d B = %.3f C = %.3f SNR = %.3f', Tb_target, B, C, SNR_cand);
fprintf('%s', str);
fprintf('\n');

% Generate modeled Rx signals
l_Pa = 10*log10(sum(A)+1); 
Pa = 10^(l_Pa / 10);
s_a = s_full;
% s_a(p.N_guard_len+p.N_sync_len:end) = s_a(p.N_guard_len+p.N_sync_len:end) ./ sqrt(Pa); %Attenuated Tx signal by distance
s_a = s_a ./ sqrt(Pa);
n_full = n;
while (length(n_full) < length(s_a))
 n_full = [n_full; n];
end
rcv_sig = (s_a + n_full(1:length(s_a)) + (20*rand(length(s_a),1) - 10));
%Add dummy silences 
rcv_sig = [n(1:p.N_guard_len,1); rcv_sig; n(1:p.N_guard_len,1)];
rcv_sig = rcv_sig  ./ 32767;

fname_rcv = fname_snd;
fname_rcv = strtok(fname_rcv,'.');
fname_rcv = [rcv_path,'/',fname_rcv,'_dist',num2str(p.dist_scale*l,'%.1f'),'.wav'];
disp(fname_rcv);
wavwrite(rcv_sig, p.fs, fname_rcv);

end
