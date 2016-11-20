% /*********************************************
% 	GIST pop-denoiser Module  (Ver.0.4m)
% 	Target: AF-noise, Pop-noise
% 	2013-08-26
% 	Human Media Communication & Processing(HuCom) Lab.
% 	
% 	Change Log.
% 				- (V0.1m) Initial release
%               - (v0.2m) Fixed perfect reconstruction condition when
%                         bypassing denoise operation
%               - (v0.2m) Applied BPF signals for both detection and
%                         reconstrudction
%               - (v0.3m) Revised reconstruction algorithm to use only
%                         previous frame
%               - (v0.4m) Tuned static tonal remover considering 300pps AF
%                         noise
%               - (v0.4m) Changed filter types from 5-th order FIR into
%                         optimized IIR
% ***********************************************/



function DRES_denoiser(prm)

%% Initialization
load('table/iir_table.mat');
input = prm.input;
output = prm.output;
FS = prm.FS;
SIZE_FRAME = prm.SIZE_FRAME;
SIZE_FSHIFT = prm.SIZE_FSHIFT;
LPorder = prm.LPorder;
BP_cutL = prm.BP_cutL;
BP_cutGAP = prm.BP_cutGAP;
BP_cutH = prm.BP_cutH;
MEAN_thr = prm.MEAN_thr;
MED_len = prm.MED_len;
MED_thr_L = prm.MED_thr_L;
MED_thr_R = prm.MED_thr_R;
TONE_len = prm.TONE_len;


[s]=wavread(input);
s=floor(s*32767);

size_total = size(s,1);
s_out = zeros(size(s));
s_frame_1d = zeros(SIZE_FRAME,1);
s_frame_bp_1d = zeros(SIZE_FRAME,1);
s_frame_bp_2d = ones(SIZE_FRAME,1);
dres_bp_1d = zeros(SIZE_FRAME,1);
detect_mask_1d = ones(SIZE_FSHIFT,1);
notch_cnt = 0;

%delay buffer for IIR filtering
dO_BP_s_frame = zeros(size(SOS_BP,1),2);
dO_BP_dres = zeros(size(SOS_BP,1),2);
dO_BP_res_noise = zeros(size(SOS_BP,1),2);
dO_BS_s_frame_1d = zeros(size(SOS_BS,1),2);
dO_notch1_s_proc = zeros(size(SOS_notch1,1),2);


size_crnt = 1;
i = 1;
w = hanning(SIZE_FRAME);
n_flag_prv = 0;
while size_crnt < size_total- SIZE_FRAME
    %% Framewise windowing
    s_frame =  w.* s(size_crnt:size_crnt-1 + SIZE_FRAME);
    
%     n = 5; %BPF tap
%     Wn = [BP_cutL BP_cutH]; %HPF cut-off
%     s_frame_bp=train_bpf(s_frame, FS, n, Wn); 
    dI_BP_s_frame = dO_BP_s_frame;
    [s_frame_bp,dO_BP_s_frame] = iir_cas(s_frame,SOS_BP,G_BP,prm.SIZE_FRAME,dI_BP_s_frame);
    
    %% Residual extraction
    LP = lpc(s_frame,LPorder);
    LP_1d = lpc(s_frame_1d,LPorder);
    res = filter(LP,1,s_frame);
    mod_res = filter(LP_1d,1,s_frame);
    
    if i > 1 %1d buffer is filled
        dres = res - mod_res;
    else     %1d buffer is empty
        dres = res;
    end
    
        
    
%     %% Detection signal generation
%     n = 5; %BPF tap
%     Wn = [BP_cutL BP_cutH]; %HPF cut-off
%     dres=train_bpf(dres, FS, n, Wn);
    dI_BP_dres = dO_BP_dres;
    [dres_bp,dO_BP_dres] = iir_cas(dres,SOS_BP,G_BP,prm.SIZE_FRAME,dI_BP_dres);
    
    dres_bp = abs(dres_bp);
    dres_bp_mid = dres_bp_1d(SIZE_FSHIFT+1:SIZE_FRAME) + dres_bp(1:SIZE_FSHIFT); 
    
    
    %% Noise area detection using mean thresholding
    detect_mask = ones(SIZE_FSHIFT,1);
    max_val = max(dres_bp_mid);
    if max_val > mean(dres_bp_mid) * MEAN_thr;
       %% Get adaptive pop_len_L and pop_len_L (check noise duration)
        pop_len_L = 0;
        pop_len_R = 0;
        for k = 1:SIZE_FSHIFT-MED_len+1
            med_crnt = median(dres_bp_mid(k:k+MED_len-1));
            
            if (pop_len_L == 0) && (max_val < med_crnt * MED_thr_L)
                pop_len_L = k;
            end
            
            if (k > pop_len_L) && (max_val > med_crnt * MED_thr_R)
                pop_len_R = k;
            end
        end
        
        if pop_len_R == 0
            pop_len_R = SIZE_FSHIFT;
        end
        
        if pop_len_L == 0
            pop_len_L = 1;
        end
        
        if pop_len_R < pop_len_L
            pop_len_R = SIZE_FSHIFT;
        end
        
        detect_mask(pop_len_L:pop_len_R) = zeros(pop_len_R-pop_len_L+1,1);
    end
    
   %% Noise region reconstruction
    
    %% Residual reconstruction
    n_flag = 0;
    n_flag_L = 0;
    n_flag_R = 0;
    if prod(detect_mask_1d) == 0
%         recon_region_L = find(~detect_mask_1d);
% 
%         recon_L_L = min(recon_region_L);
%         recon_L_R = max(recon_region_L);
%         
%         dur_recon_L = recon_L_R - recon_L_L + 1;
%         
          n_flag = 1;
%          n_flag_L = 1;
    end
    
    if prod(detect_mask) == 0
%         recon_region_R = find(~detect_mask);
% 
%         recon_R_L = min(recon_region_R)+SIZE_FSHIFT;
%         recon_R_R = max(recon_region_R)+SIZE_FSHIFT;
%       
%         dur_recon_R = recon_R_R - recon_R_L + 1;
% 
           n_flag = 1;
%          n_flag_R = 1;
    end
    
    if prm.BYPASS == 1
        n_flag = 0;
    end
    
    if prm.BYPASS_INV == 1
        n_flag = 1;
    end
    
    %% Frame Re-synthesis
    if n_flag == 0 %|| ~(sum(abs(s_frame_rect_2d)) > 0)
        s_proc = s_frame_1d;
        s_frame_bp_2d = s_frame_bp_1d;
    else
        %LP-based waveform resynthesis (source-filter)
        LPorder_rec_L = prm.LPorder_recon;

        [LP_rec_b] = lpc(s_frame_bp_2d(SIZE_FSHIFT+1:SIZE_FRAME), LPorder_rec_L);

        
        %Residual generation
        %====TABLE====
        res_noise = (2*rand(SIZE_FRAME,1)-1); 
        dI_BP_res_noise = dO_BP_res_noise;
        [res_noise,dO_BP_res_noise] = iir_cas(res_noise,SOS_BP,G_BP,prm.SIZE_FRAME,dI_BP_res_noise);
        %====TABLE====
        res_bp = filter(LP_rec_b,1,s_frame_bp_2d(SIZE_FSHIFT+1:SIZE_FRAME));
        res_pow = mean(abs(res_bp));
        res_noise = res_noise .* res_pow; %residual power matching
        
        
%         n = 5; %BPF tap
%         Wn = [BP_cutL BP_cutH]; %HPF cut-off
%         res_noise=train_bpf(res_noise, FS, n, Wn);    

        %Synthesize reconstruction frame
        s_proc_bp = filter(1, LP_rec_b, res_noise);
        s_proc_bp = s_proc_bp * prm.RECON_POWER;
        
        
%         n = 5; %BPF tap
%         Wn = [BP_cutL BP_cutH]; %HPF cut-off
%         s_proc_bs=train_bsf(s_frame_1d, FS, n, Wn);
        dI_BS_s_frame_1d = dO_BS_s_frame_1d;
        [s_proc_bs,dO_BS_s_frame_1d] = iir_cas(s_frame_1d,SOS_BS,G_BS,prm.SIZE_FRAME,dI_BS_s_frame_1d);
        
        s_proc = s_proc_bp + s_proc_bs;
        
        if prm.DETECT == 1
            s_proc = s_proc_bs;
        end
        
        s_frame_bp_2d = s_proc_bp;
    end
    
    if prm.TONE_REMOVE == 1
        %% remove tonal noise using notch filter
        if n_flag_prv == 1 || (notch_cnt > 0 && notch_cnt <= TONE_len)
            %s_proc = train_bsf(s_proc, FS, 3, [5900 6300]); 
            dI_notch1_s_proc = dO_notch1_s_proc;
            [s_proc,dO_notch1_s_proc] = iir_cas(s_proc,SOS_notch1,G_notch1,prm.SIZE_FRAME,dI_notch1_s_proc);
            
            notch_cnt = notch_cnt+1;
        end

        if notch_cnt == TONE_len
            notch_cnt = 0;
        end
    end
    
    %% DRES write option
    if prm.DRES_OUT == 1
        s_proc = dres_bp;
    end
    
    s_out(size_crnt : size_crnt-1 + SIZE_FRAME) = s_out(size_crnt : size_crnt-1 + SIZE_FRAME) + s_proc;
    
    
    %% Proccessed signal buffering
    s_frame_1d = s_frame;
    s_frame_bp_1d = s_frame_bp;
    dres_bp_1d = dres_bp;
    detect_mask_1d = detect_mask;
    n_flag_prv = n_flag;
    
    size_crnt = size_crnt + SIZE_FSHIFT;
    i = i+1;
end

wavwrite(s_out./32767, FS, 16, output);

% plot(s);
% hold on;
% plot(s_out(SIZE_FSHIFT:length(s_out))*100,'r');

end