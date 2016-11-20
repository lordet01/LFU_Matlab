% clear;
% clc;

mkdir('input');
mkdir('output');
prm.input = 'wav/AC_smooth_c_MP.wav';      
prm.output = 'wav/AC_smooth_c_MP_detect.wav';
prm.FS = 44100;     
prm.SIZE_FRAME = 1024;
prm.SIZE_FSHIFT = 1024;
prm.LPorder = 10;  
prm.MEAN_thr = 50;

prm.DRES_OUT =0;
prm.DETECT = 1;

addpath('src');


input = prm.input;
output = prm.output;

[s]=wavread(input);
s=floor(s*32767);

[peak_pos, flag, s_out] = DRES_detect(s, prm);

flag = [flag(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)];
wavwrite([s_out(prm.SIZE_FRAME * 2+1:end); zeros(prm.SIZE_FRAME * 2,1)]./32767, prm.FS, 16, output);

plot(s);
hold on;
plot(flag,'ro');
% peak_buff = zeros(size(flag));
% peak_buff(peak_pos - prm.SIZE_FRAME * 4) = max(s);
% plot(peak_buff, 'go');
 
