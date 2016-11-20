clc; clear; close('all')

load('BPF_44100_18050_22050');
BPF_tab = length(BPF);
fs = 44100;
t1 = 0.2; %400ms
fl = 18050;
fh = 22050;
ws = 1.0; % > 1 thinner, <1 : thicker
rr = 0.9; %Rectification ratio

%quadratic chirp
t = -t1:1/fs:t1;
y1 = chirp(t,fl,t1,fh,'quadratic',[],'concave');
y2 = chirp(t,fl,t1,fh,'quadratic',[],'convex');

t_h = floor(length(t) / 2);
t_win = floor(t_h*(1-rr));
w_tail = hann(t_win);
w_tail = w_tail .^ ws;
w = ones(t_h,1);
w(1:floor(t_win/2)) = w_tail(1:floor(t_win/2));
w(end-floor(t_win/2)+1:end) = w_tail(1+floor(t_win/2):2*floor(t_win/2));
y1 = y1(1:t_h) .* w';
y2 = y2(t_h+2:end) .* w';

p = [y1 , y2];  

% % % linear chirp
% % t = 0:1/fs:t1;
% % y2 = chirp(t,fl,t1,fh,'linear');
% % y1 = chirp(t,fh,t1,fl,'linear');
% % 
% % t_h = floor(length(t));
% % w = hann(t_h);
% % w = w .^ ws;
% % y1 = y1(1:t_h) .* w';
% % y2 = y2(1:t_h) .* w';
% % 
% % p = [y1 , y2];  


% p = conv(p,BPF);

% spectrogram(y,2048,2044,2048,fs,'yaxis')
wavwrite(p, fs, 'preamble.wav');


%% Multipath profile (16-MOBICOM, Soundlly)
Td = [0, 0.0002, 0.0005, 0.01, 0.025, 0.03]; %s
Att = [0, -3.78, -2.13, -2.40, -4, -5.89]; %dB
Td_s = ceil(fs.*Td) + 1;
Att_g = 1 ./ (10.^(-1 * Att / 20));
H = zeros(Td_s(end), 1);
H(Td_s) = Att_g;
% H = 1.0; %Turn off Multipath effect

%Filter the multipath
y = conv(H, p);
y = y / max(y) * 0.2;


%% Add background noise
SNR = 0; %Subband SNR %max: 60
n = wavread('noise_218.wav');
pre_start = 7000;
n_sub = n(pre_start+1:pre_start+length(y));
g_n = sqrt(mean(conv(y,BPF).^2) / mean(conv(n_sub,BPF).^2) / (10^(SNR/10)));

% n = n * g_n; n_sub = n_sub * g_n;
n(pre_start+1:pre_start+length(y)) = n_sub + y' / g_n;
y = n;
wavwrite(y, fs, 'preamble_sim_c_MP.wav');


%% Ideal Rake reciever test 

%Multipath compensation test (what if we know lags of multipath?)
H_nz_I=find(H>0);
MP_min = 50;
y_sum = y;
for k = 2:length(H_nz_I)
    if H_nz_I(k) > MP_min
%      H_g = H(H_nz(k));
     H_g = 1;
     y_sum = y_sum + [H_g*y(H_nz_I(k):end); zeros(H_nz_I(k) - 1,1)];
    end
end
wavwrite(y_sum, fs, 'preamble_sim_c_MP_Irake.wav');

% % y = y_sum; 

%% Correlation Test (Matched Filtering)

beta = 4; %Power factor for correlated result

y_BPF = conv(y,BPF);
AC = conv(y_BPF, fliplr(p));
AC = AC.^beta;

%smoothing % scaling
lag = 0;
if lag > 0
    AC_smooth = tsmovavg(AC,'s',lag,1);
    AC_smooth = AC_smooth(lag+1:end);
else
    AC_smooth = AC;
end

AC_smooth = AC_smooth + 0.00001; %prevent zeros
wavwrite(AC_smooth ./ max(AC_smooth),fs,'AC_smooth_c_MP.wav');
% plot(AC_smooth,'o');

%% Get First peak position
prm.SIZE_FRAME = floor(t1*fs);
prm.SIZE_FSHIFT = floor(t1*fs);
prm.LPorder = 10;  
prm.MEAN_thr = 500;
prm.DRES_OUT =0;
prm.DETECT = 1;

[peak_pos, flag, s_out] = DRES_detect(AC_smooth, prm);

flag = [flag(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)];
% wavwrite([s_out(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)]./32767, prm.FS, 16, output);

figure();
plot(AC_smooth);
hold on;
plot(flag,'ro');
hold off;

% sync accuracy
peak_pos = find(flag>0.001);
peak_pos = peak_pos(1);
pre_start_est = peak_pos - floor((BPF_tab -1)*0.5) - length(p);
acc = (pre_start_est - pre_start) / fs;
fprintf('%.3f ms\n',acc * 1000);


%% Rake reciever  

%search fingers 
Td_pool = 0.1; %100ms
Td_pool_s = ceil(fs.*Td_pool) + 1;
Td_gap = 0.005; %10ms
Td_gap_s = ceil(fs.*Td_gap) + 1;
AC_sub = AC_smooth(peak_pos: peak_pos+Td_pool_s);

prm.SIZE_FRAME = Td_gap_s;
prm.SIZE_FSHIFT = Td_gap_s;
prm.MEAN_thr = 1;
[~, flag, ~] = DRES_detect(AC_sub, prm);

flag = [flag(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)];
% wavwrite([s_out(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)]./32767, prm.FS, 16, output);

figure();
plot(AC_sub);
hold on;
plot(flag,'ro');
hold off;


% sum-up multipath signals
H_nz=find(flag>flag(1)*0.05);
H_nz = H_nz +1; % Need to search optimum value (around -3~+3)
MP_min = 50;
H_est = AC_sub ./ max(AC_sub);
H_est = H_est .^ (1/(beta*2));
H_g = H_est(H_nz);
H_g = 1;
y_sum = y;
for k = 2:length(H_nz)
    if H_nz(k) > MP_min
% %      y_finger= y;
% %      y_finger(1:pre_start_est+H_nz(k)-1) = zeros(pre_start_est+H_nz(k)-1, 1);
% %      y_finger = [H_g(k)*y_finger(H_nz(k):end); zeros(H_nz(k) - 1,1)];
% %      
% %      y_sum = y_sum + y_finger;
     y_sum = y_sum + [H_g*y(H_nz(k):end); zeros(H_nz(k) - 1,1)];
    end
end

wavwrite(y_sum, fs, 'preamble_sim_c_MP_rake.wav');



%% Correlation Test (Matched Filtering)

beta = 4; %Power factor for correlated result

y_BPF = conv(y_sum,BPF);
AC = conv(y_BPF, fliplr(p));
AC = AC.^beta;

%smoothing % scaling
lag = 0;
if lag > 0
    AC_smooth = tsmovavg(AC,'s',lag,1);
    AC_smooth = AC_smooth(lag+1:end);
else
    AC_smooth = AC;
end

AC_smooth = AC_smooth + 0.00001; %prevent zeros
wavwrite(AC_smooth ./ max(AC_smooth),fs,'AC_smooth_c_MP.wav');
% plot(AC_smooth,'o');

%% Get First peak position
prm.SIZE_FRAME = floor(t1*fs);
prm.SIZE_FSHIFT = floor(t1*fs);
prm.LPorder = 10;  
prm.MEAN_thr = 500;
prm.DRES_OUT =0;
prm.DETECT = 1;

[peak_pos, flag, s_out] = DRES_detect(AC_smooth, prm);

flag = [flag(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)];
% wavwrite([s_out(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)]./32767, prm.FS, 16, output);

figure();
plot(AC_smooth);
hold on;
plot(flag,'ro');
hold off;

% sync accuracy
peak_pos = find(flag>0.001);
peak_pos = peak_pos(1);
peak_pos = peak_pos - floor((BPF_tab -1)*0.5);
pre_start_est = peak_pos - length(p);
acc = (pre_start_est - pre_start) / fs;
fprintf('%.3f ms\n',acc * 1000);


% % % % % 
% % % % % 
% % % % % % % Get results
% % % % % % y = AC;
% % % % % % [signals,avg,dev] = threshold_zscore(y,lag,threshold,influence);
% % % % % % 
% % % % % % figure; subplot(2,1,1); hold on;
% % % % % % x = 1:length(y); ix = lag+1:length(y);
% % % % % % area(x(ix),avg(ix)+threshold*dev(ix),'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
% % % % % % area(x(ix),avg(ix)-threshold*dev(ix),'FaceColor',[1 1 1],'EdgeColor','none');
% % % % % % plot(x(ix),avg(ix),'LineWidth',1,'Color','cyan','LineWidth',1.5);
% % % % % % plot(x(ix),avg(ix)+threshold*dev(ix),'LineWidth',1,'Color','green','LineWidth',1.5);
% % % % % % plot(x(ix),avg(ix)-threshold*dev(ix),'LineWidth',1,'Color','green','LineWidth',1.5);
% % % % % % plot(1:length(y),y,'b');
% % % % % % subplot(2,1,2);
% % % % % % stairs(signals,'r','LineWidth',1.5); ylim([-1.5 1.5]);
% % % % % 




