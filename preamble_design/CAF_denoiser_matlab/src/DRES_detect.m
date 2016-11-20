% /*********************************************
% 	GIST pop detector
% 	2016-11-16
% 	Human Media Communication & Processing(HuCom) Lab.
%
% 	Change Log.
% 				- (V0.1m) Initial release (modified from GIST pop denoiser
% 				  in Samsung Project 2013)
% ***********************************************/



function [peak_pos, flag_out, s_out]=DRES_detect(s, prm)

%% Initialization
SIZE_FRAME = prm.SIZE_FRAME;
SIZE_FSHIFT = prm.SIZE_FSHIFT;
LPorder = prm.LPorder;
MEAN_thr = prm.MEAN_thr;


size_total = size(s,1);
s_out = zeros(size(s));
flag_out = zeros(size(s));

%Temporal buffers for operation
s_frame_1d = zeros(SIZE_FRAME,1) + 2;
s_frame_LP_1d = zeros(SIZE_FRAME*2,1) + 2;
dres_bp_1d = zeros(SIZE_FRAME*2,1) + 0.000001;
% n_flag_prv = 0;



frame_idx = 1;
cnt = 1;


size_crnt = 1;
i = 1;
w = hanning(SIZE_FRAME*2);
% w = ones(SIZE_FRAME*2,1);
while size_crnt < size_total- SIZE_FRAME*2
    %% Framewise windowing
    s_frame = s(size_crnt:size_crnt-1 + SIZE_FRAME);
    s_frame_LP =  w.* [s_frame_1d ; s_frame];

    %% Residual extraction
    LP = lpc(s_frame_LP,LPorder);
    LP_1d = lpc(s_frame_LP_1d,LPorder);
    res = filter(LP,1,s_frame_LP);
    mod_res = filter(LP_1d,1,s_frame_LP);

    if i > 1 %1d buffer is filled
        dres = res - mod_res;
    else     %1d buffer is empty
        dres = res;
    end
    
    dres_bp = dres;
    dres_bp = abs(dres_bp);
    dres_bp_mid = dres_bp_1d(SIZE_FSHIFT+1:SIZE_FRAME*2) + dres_bp(1:SIZE_FRAME*2 - SIZE_FSHIFT);

    %% Peak area detection by mean-max thresholding
    max_val = max(dres_bp_mid);
    mean_val = mean(dres_bp_mid);
    cnt = cnt+1;

    if max_val > mean_val * MEAN_thr
        n_flag = 1;
    else
        n_flag = 0;
    end
    

    if n_flag == 1 %|| n_flag_prv == 1
        n_flag_final = ones(SIZE_FRAME*2,1);
        s_proc = s_frame_LP_1d .* 0;
        [peak_val,peak_pos] = max(s_frame_LP_1d);
        n_flag_final(peak_pos) = peak_val;
    else
        n_flag_final = zeros(SIZE_FRAME*2,1);
        s_proc = s_frame_LP_1d;
        peak_pos = 0;
    end

    %% DRES write option
    if prm.DRES_OUT == 1
        s_proc = dres_bp;
    end
    
    s_out(size_crnt : size_crnt-1 + SIZE_FRAME*2) = s_out(size_crnt : size_crnt-1 + SIZE_FRAME*2) + s_proc;
    flag_out(size_crnt : size_crnt-1 + SIZE_FRAME*2) = flag_out(size_crnt : size_crnt-1 + SIZE_FRAME*2) + n_flag_final;

    %% Proccessed signal buffering
    s_frame_1d = s_frame;
    s_frame_LP_1d = s_frame_LP;
    dres_bp_1d = dres_bp; 
    
%     n_flag_prv = n_flag;

    size_crnt = size_crnt + SIZE_FSHIFT;
    i = i+1;
    frame_idx = frame_idx+1;
end


% plot(s);
% hold on;
% plot(s_out(SIZE_FSHIFT:length(s_out))*100,'r');

end