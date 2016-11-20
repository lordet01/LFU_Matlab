% /*********************************************
% 	GIST pop-denoiser Module  (Ver.0.5m)
% 	Target: AF-noise, Pop-noise
% 	2013-08-27
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
%               - (v0.5m) Revised frame reconstruction algorithm based on
%                         MDCT analysis
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
s_frame_rect_1d = zeros(SIZE_FRAME,1);
S_in = zeros(SIZE_FRAME*0.5,1);
S_in_1d = zeros(SIZE_FRAME*0.5,1);

dres_bp_1d = zeros(SIZE_FRAME,1);
detect_mask_1d = ones(SIZE_FSHIFT,1);
notch_cnt = 0;

%delay buffer for IIR filtering
% dO_BP_s_frame = zeros(size(SOS_BP,1),2);
dO_BP_dres = zeros(size(SOS_BP,1),2);
% dO_BP_res_noise = zeros(size(SOS_BP,1),2);
% dO_BS_s_frame_1d = zeros(size(SOS_BS,1),2);
dO_notch1_s_proc = zeros(size(SOS_notch1,1),2);


size_crnt = 1;
i = 1;
w = hanning(SIZE_FRAME);
n_flag_prv = 0;
while size_crnt < size_total- SIZE_FRAME
    %% Framewise windowing
    s_frame =  w.* s(size_crnt:size_crnt-1 + SIZE_FRAME);
    s_frame_rect = s(size_crnt:size_crnt-1 + SIZE_FRAME);
    

%     dI_BP_s_frame = dO_BP_s_frame;
%     [s_frame_bp,dO_BP_s_frame] = iir_cas(s_frame,SOS_BP,G_BP,prm.SIZE_FRAME,dI_BP_s_frame);
    
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
    S_in = mdct(s_frame_rect_1d);
    if n_flag == 1 %|| ~(sum(abs(s_frame_rect_2d)) > 0)
        %Magnitude reconstruction
        binL = floor((BP_cutL / FS) * SIZE_FRAME);
        binH = floor((BP_cutH / FS) * SIZE_FRAME);
        S_in(binL:binH) = S_in_1d(binL:binH) .* sign(2.*rand(size(binL:binH))'-1);
        
        if prm.DETECT == 1
            S_in(binL:binH) = ones((binL:binH),1)*10000;
        end
    end
    
    if prm.TONE_REMOVE == 1
        %% remove tonal noise using notch filter
        if n_flag_prv == 1 || (notch_cnt > 0 && notch_cnt <= TONE_len)
            %dI_notch1_s_proc = dO_notch1_s_proc;
            %[s_proc,dO_notch1_s_proc] = iir_cas(s_proc,SOS_notch1,G_notch1,prm.SIZE_FRAME,dI_notch1_s_proc);
            binL_notch1 = floor((5900 / FS) * SIZE_FRAME);
            binH_notch1 = floor((6300 / FS) * SIZE_FRAME);
            S_in(binL_notch1:binH_notch1) =  S_in_1d(binL_notch1:binH_notch1)  .* sign(2.*rand(size(binL_notch1:binH_notch1))'-1);
            
            notch_cnt = notch_cnt+1;
        end

        if notch_cnt == TONE_len
            notch_cnt = 0;
        end
    end
    
    s_proc = imdct(S_in);

    
    %% DRES write option
    if prm.DRES_OUT == 1
        s_proc = dres_bp;
    end
    
    s_out(size_crnt : size_crnt-1 + SIZE_FRAME) = s_out(size_crnt : size_crnt-1 + SIZE_FRAME) + s_proc;
    
    
    %% Proccessed signal buffering
    S_in_1d = S_in;
    s_frame_1d = s_frame;
    s_frame_rect_1d = s_frame_rect;
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