

function [A_prm]=atten_opt(Tx_file, real_path, Tb, fc, p)

%%Parameters
fs = p.fs;
g = p.g; %Spreading factor
% As = p.As; %Unit-normalizing constant for Attenuation factor

%Load Tx files
s_Tx = wavread(Tx_file);
s_Tx = s_Tx .* 32767;
s_Tx = s_Tx(p.comp_len+1:2*p.comp_len);
K = length(s_Tx);
K21 = K/2+1;
% fc_bin = floor(fc / fs * 2 * K21);

fc_L = fc - fs/Tb;
fc_H = fc + fs/Tb;
fc_bin = fc / fs * 2 * K21;
%Calculate band of interest
B_init_f = (fc_L : fc_H);
B_init_bin = floor(B_init_f ./ fs .* 2 .* K21);


%% Get Eb/N0 for given l and SNRs
%Get PSD of Tx
S_Tx = abs(fft(s_Tx, K));
S_Tx = sum(S_Tx(B_init_bin).^2);


%% Reference attenuation 
ref_list = dir(real_path);
ydata = zeros(1, length(ref_list)-2);
for i = 3:length(ref_list)
    
    %Calculate real attenuation
    l = i-2;
    s_real = wavread([real_path,'/',ref_list(i).name]);
%     disp(ref_list(i).name);
    s_real = s_real .* 32767;
    s_real = s_real(p.comp_len+1:2*p.comp_len);
    S_real = abs(fft(s_real, K));
    S_real = sum(S_real(B_init_bin).^2); %./ (K*fs);
%     disp(S_real);
    ydata(l) = S_Tx ./ S_real; %Attenuation ratio of real data
    
end
    %Attenuation function
%     B_init_f = B_init_bin / (K21) * fs * pi;
%     A = As * l^g * exp(l * (alpha_0 * B_init_f .^ eta));
    
%find optimum parameters by least-square curve fitting
A0 = [1, 0.001, 0.9];
A_LB = [0, 0, 0]; %Lower bound
R = p.R_dist;
% xdata = [1/R:1/R:l]';
xdata = [1:1:l*R]';
ydata = sort(ydata(:),'ascend');
ydata = interp(ydata, R, 2);
funA = @(x,xdata)(x(1) * xdata.^g) .* exp(xdata .* (x(2) * fc_bin .^ x(3)));
A_prm = lsqcurvefit(funA,A0,xdata,ydata, A_LB);
% As = x(1);
% alpha_0 = x(2);
% eta = x(3);

% display fitting results
times = linspace(xdata(1),xdata(end));
plot(xdata,ydata,'ko',times,funA(A_prm,times),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')





