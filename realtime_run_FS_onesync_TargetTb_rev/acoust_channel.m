function y = acoust_channel(x, d, fs)

load('channel/BPF_44100_18050_22050'); 
BPF=BPF'; BPF_tab = length(BPF);

LOS_COND = 2; %0: IDEAL, 1: LOS, 2: NLOS
START_POS = 60000; %Start position of comm. signal (consider noise signal is always longer!)
SNR = -1; %dB Scale

% %% Multipath profile (16-MOBICOM, Soundlly)
if LOS_COND == 1
    %LOS
    Td = [0, 0.0002, 0.0005, 0.01, 0.025, 0.03]; %s
    Att = [0.0, -3.78, -2.13, -2.40, -4, -5.89]; %dB
    Td_s = ceil(fs.*Td) + 1;
    Att_g = 1 ./ (10.^(-1 * Att / 20));
    H = zeros(Td_s(end), 1);
    H(Td_s) = Att_g;
    
elseif LOS_COND == 2
    %NLOS
    Td = [0, 0.0002, 0.0005, 0.01, 0.025, 0.03]; %s
    Att = [-2.11, -3.78, -2.13, -2.40, -4.0, -5.89]; %dB
    Td_s = ceil(fs.*Td) + 1;
    Att_g = 1 ./ (10.^(-1 * Att / 20));
    H = zeros(Td_s(end), 1);
    H(Td_s) = Att_g;
else
    H = 1;
end

% % %Filter the multipath (Consider received preamble)
x_hat = conv(H, x);

%% Add background noise
pre_start = START_POS;

d_sub = d(pre_start+1:pre_start+length(x_hat));
g_n = sqrt(mean(conv(x_hat,BPF).^2) / mean(conv(d_sub,BPF).^2) / (10^(SNR/10)));
y = d;
y(pre_start+1:pre_start+length(x_hat)) = d_sub + x_hat / g_n;


end