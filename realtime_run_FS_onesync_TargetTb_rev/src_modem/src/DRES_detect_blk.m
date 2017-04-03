% /*********************************************
% 	DRES preamble detector (block-wise)
% 	2017-04-04
% 	Human Media Communication & Processing(HuCom) Lab.
%
% 	Change Log.
% 				- (V1.0m) Initial release (modified from GIST pop denoiser
% 				  in Samsung Project 2013)
% ***********************************************/



function [peak_pos, g]=DRES_detect_blk(s_blk, g, prm)

%% Initialization
SIZE_FRAME = prm.SIZE_FRAME;
SIZE_FSHIFT = prm.SIZE_FSHIFT;
LPorder = prm.LPorder;
MEAN_thr = prm.MEAN_thr;


%% Temporal buffer updates (input)
s_frame_1d = g.s_frame_1d;
s_frame_LP_1d = g.s_frame_LP_1d;
peak_val_1d = g.peak_val_1d;
dres_bp_1d = g.dres_bp_1d;
cnt = g.cnt;
frame_idx = g.frame_idx;
w = g.w;
peak_val_sel = g.peak_val_sel;

    
    %% Framewise windowing
    %Initialize 1d buffers
    if frame_idx == 1
         s_frame_1d = s_blk;
    end
    
    s_frame_LP =  w.* [s_frame_1d ; s_blk];
    % Initialize 1d buffers
    if frame_idx == 1
         s_frame_LP_1d = s_frame_LP;
    end
    
    %% Residual extraction
    LP = lpc(s_frame_LP,LPorder);
    LP_1d = lpc(s_frame_LP_1d,LPorder);
    res = filter(LP,1,s_frame_LP);
    mod_res = filter(LP_1d,1,s_frame_LP);

    if frame_idx > 1 %1d buffer is filled
        dres = res - mod_res;
    else     %1d buffer is empty
        dres = res;
    end
    
    dres_bp = dres;
    dres_bp = abs(dres_bp);
    
    if frame_idx == 1
         dres_bp_1d = dres_bp;
    end
    
    dres_bp_mid = dres_bp_1d(SIZE_FSHIFT+1:SIZE_FRAME*2) + dres_bp(1:SIZE_FRAME*2 - SIZE_FSHIFT);

    %% Peak area detection by mean-max thresholding
    max_val = max(dres_bp_mid);
    mean_val = mean(dres_bp_mid);
    cnt = cnt+1;
    peak_pos = 0; %default output
    
    if max_val > mean_val * MEAN_thr && mean_val > 0.000001
        [peak_val,peak_pos_cand] = max(dres_bp_1d);
        if peak_val_1d > 0.001 && peak_val > peak_val_1d * MEAN_thr && peak_val > peak_val_sel
            peak_pos = peak_pos_cand;
            peak_val_sel = peak_val;
        end
        peak_val_1d = peak_val;
    end
    
    
%% Temporal buffer updates (output)
g.s_frame_1d = s_blk;
g.s_frame_LP_1d = s_frame_LP;
g.peak_val_1d = peak_val_1d;
g.dres_bp_1d = dres_bp;
g.frame_idx = frame_idx;
g.cnt = cnt;
g.peak_val_sel = peak_val_sel;
    
    
end
