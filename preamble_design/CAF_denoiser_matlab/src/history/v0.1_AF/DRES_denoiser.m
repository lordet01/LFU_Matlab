% /*********************************************
% 	GIST pop-denoiser Module  (Ver.0.1m)
% 	Target: AF-noise, Pop-noise
% 	2013-06-14
% 	Human Media Communication & Processing(HuCom) Lab.
% 	
% 	Change Log.
% 				- (V0.1m) Initial release
% ***********************************************/



function DRES_denoiser(prm)

%% Initialization
input = prm.input;
output = prm.output;
fs = prm.fs;
size_frame = prm.size_frame;
size_fshift = prm.size_fshift;
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
s_frame_1d = zeros(size_frame,1);
s_frame_2d = zeros(size_frame,1);
s_frame_rect_1d = zeros(size_frame,1);
s_frame_rect_2d = zeros(size_frame,1);
res_1d = zeros(size_frame,1);
dres_bp_1d = zeros(size_frame,1);
detect_mask_1d = ones(size_fshift,1);
LP_2d = zeros(1,LPorder+1); LP_2d(1) = 1;
notch_cnt = 0;

size_crnt = 1;
i = 1;
w = hanning(size_frame);
n_flag_prv = 0;
while size_crnt < size_total- size_frame
    %% Framewise windowing
    s_frame =  w.* s(size_crnt:size_crnt-1 + size_frame);
    s_frame_rect = s(size_crnt:size_crnt-1 + size_frame);
    
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
    
        
    
    %% Detection signal generation
    n = 5; %BPF tap
    Wn = [BP_cutL BP_cutH]; %HPF cut-off
    dres_bp=train_bpf(dres, fs, n, Wn);
    
    dres_bp = abs(dres_bp);
    
    dres_bp_mid = dres_bp_1d(size_fshift+1:size_frame) + dres_bp(1:size_fshift); 
    
    
    %% Noise area detection using mean thresholding
    detect_mask = ones(size_fshift,1);
    max_val = max(dres_bp_mid);
    if max_val > mean(dres_bp_mid) * MEAN_thr;
       %% Get adaptive pop_len_L and pop_len_L (check noise duration)
        pop_len_L = 0;
        pop_len_R = 0;
        for k = 1:size_fshift-MED_len+1
            med_crnt = median(dres_bp_mid(k:k+MED_len-1));
            
            if (pop_len_L == 0) && (max_val < med_crnt * MED_thr_L)
                pop_len_L = k;
            end
            
            if (k > pop_len_L) && (max_val > med_crnt * MED_thr_R)
                pop_len_R = k;
            end
        end
        
        if pop_len_R == 0
            pop_len_R = size_fshift;
        end
        
        if pop_len_L == 0
            pop_len_L = 1;
        end
        
        if pop_len_R < pop_len_L
            pop_len_R = size_fshift;
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
%         recon_R_L = min(recon_region_R)+size_fshift;
%         recon_R_R = max(recon_region_R)+size_fshift;
%       
%         dur_recon_R = recon_R_R - recon_R_L + 1;
% 
           n_flag = 1;
%          n_flag_R = 1;
    end
    
    
    %% Frame Re-synthesis
    if n_flag == 0 %|| ~(sum(abs(s_frame_rect_2d)) > 0)
        s_proc = s_frame_1d;
        
        s_frame_rect_2d = s_frame_rect_1d;
    else
        %LP-based waveform resynthesis (source-filter)
        LPorder_rec_L = 500;
        LPorder_rec_R = 200;
        
        [LP_rec_b, var_L] = lpc(s_frame_rect_2d(1:size_fshift), LPorder_rec_L);
        [LP_rec_f, var_R] = lpc(s_frame_rect(size_fshift+1+200:size_frame), LPorder_rec_R);

        res_noise_L = (2*rand(size_frame,1)+1);
        res_noise_R = (2*rand(size_frame,1)+1);
        w2 = hanning(2*size_frame);
        s_proc_bp_rec_L = sqrt(var_L) * filter(1, LP_rec_b, res_noise_L);
        s_proc_bp_rec_R = sqrt(var_R) * filter(1, LP_rec_f, res_noise_R);

        %s_proc_bp_rec = s_proc_bp_rec_L;
        s_proc_bp_rec_L = s_proc_bp_rec_L .* w2(size_frame+1:2*size_frame);
        s_proc_bp_rec_R = s_proc_bp_rec_R .* w2(1:size_frame);
        
        s_proc_bp_rec = (s_proc_bp_rec_L + s_proc_bp_rec_R);
        
%         if n_flag_L == 1 && n_flag_R == 1
%             s_proc_bp_rec = [s_frame_rect_1d(1:recon_L_L) ; s_proc_bp_rec(recon_L_L+1:recon_R_R) ; s_frame_rect_1d(recon_R_R+1:size_frame)];
%         elseif n_flag_L == 1 && n_flag_R == 0
%             s_proc_bp_rec = [s_frame_rect_1d(1:recon_L_L) ; s_proc_bp_rec(recon_L_L+1:recon_L_R) ; s_frame_rect_1d(recon_L_R+1:size_frame)];
%         elseif  n_flag_L == 0 && n_flag_R == 1
%             s_proc_bp_rec = [s_frame_rect_1d(1:recon_R_L) ; s_proc_bp_rec(recon_R_L+1:recon_R_R) ; s_frame_rect_1d(recon_R_R+1:size_frame)];
%         end
        
        n = 5; %BPF tap
        Wn = [BP_cutL BP_cutH]; %HPF cut-off
        s_proc_bp_rec=train_bpf(s_proc_bp_rec, fs, n, Wn);        
        s_frame_rect_2d = train_bpf(s_frame_rect_2d, fs, n, Wn);    
        s_frame_rect = train_bpf(s_frame_rect, fs, n, Wn);    
        
        mean_b = mean(abs(s_frame_rect_2d(1 : size_fshift)));
        mean_f = mean(abs(s_frame_rect(size_fshift+1 : size_frame)));
        mean_c = mean(abs(s_proc_bp_rec)) * 1.5;
        
        tilt = (mean_f - mean_b) / (size_frame-1);
        mean_curve = (1:size_frame)';
        mean_curve = (mean_curve -1) .* tilt + mean_b; 
        s_proc_bp_rec = s_proc_bp_rec .* mean_curve / mean_c;

%         %LP-based waveform resynthesis (AR model)
%         LPorder_rec = 1000;
%         [LP_rec_b, var_L] = lpc(s_frame_rect_2d,LPorder_rec);
%         [LP_rec_f, var_R] = lpc(s_frame_rect,LPorder_rec);
%         s_proc_bp_rec = zeros(size_frame,1);
%         s_proc_bp_rec_L = zeros(size_frame,1);
%         s_proc_bp_rec_R = zeros(size_frame,1);
%         s_frame_buff_L = [s_frame_rect_2d; s_frame_rect_1d; s_frame_rect];
%         s_frame_buff_R = [s_frame_rect_2d; s_frame_rect_1d; s_frame_rect];
%         
%         for j = 1 : size_frame
%             for lp = 2:LPorder_rec
%                 %Left to right reconstruction
%                 s_proc_bp_rec_L(j) = s_proc_bp_rec_L(j) -1 * LP_rec_b(lp) * s_frame_buff_L(size_frame+j-(lp-1));
%                 
%                 %Right to left reconstruction
%                 s_proc_bp_rec_R(size_frame-j+1) = s_proc_bp_rec_R(size_frame-j+1) -1 * LP_rec_f(lp) * s_frame_buff_R(2*size_frame-j+(lp-1)+1);
%             end
%             s_frame_buff_L(j+size_frame) = s_proc_bp_rec_L(j);
%             s_frame_buff_R(2*size_frame-j+1) = s_proc_bp_rec_R(size_frame-j+1);
%         end
%         s_proc_bp_rec = (s_proc_bp_rec_L + s_proc_bp_rec_R) * 0.5;

        
%         %LP-based spectrum resynthesis
%         hop_size = 128;
%         FFT_len = 512; %< size_frame
%         s_frame_buff_L = [s_frame_rect_2d; zeros(size_frame,1); s_frame_rect];
%         s_frame_buff_R = [s_frame_rect_2d; zeros(size_frame,1); s_frame_rect];
%         
%         spec_buff_len = (size_frame - size_fshift) / hop_size;
%         spec_recon_len = size_fshift / hop_size;
%         spec_L = zeros(FFT_len*0.5+1, spec_buff_len+spec_recon_len);
%         spec_R = zeros(FFT_len*0.5+1, spec_buff_len+spec_recon_len);
%         sub_win = hanning(FFT_len);
%         for fft_center = 1: spec_buff_len
%             fft_idx_L = size_fshift + hop_size * (fft_center-1);
%             fft_mag_L = abs(fft(s_frame_buff_L(fft_idx_L-FFT_len*0.5 : fft_idx_L+FFT_len*0.5-1) .* sub_win));
%             spec_L(:,fft_center) = fft_mag_L(1:FFT_len*0.5+1);
%             
%             fft_idx_R = 2*size_frame + size_fshift - hop_size * (fft_center-1);
%             fft_mag_R = abs(fft(s_frame_buff_R(fft_idx_R-FFT_len*0.5 : fft_idx_R+FFT_len*0.5-1) .* sub_win));
%             spec_R(:,fft_center) = fft_mag_R(1:FFT_len*0.5+1);
%         end
%         
%         bin_L = floor(FFT_len / fs * BP_cutL);
%         bin_H = ceil(FFT_len / fs * BP_cutH);
%         
%         LP_spec_order = spec_buff_len-1; % should be < spec_buff_len
%         LP_spec_L = zeros(LP_spec_order+1, FFT_len*0.5 + 1);
%         LP_spec_R = zeros(LP_spec_order+1, FFT_len*0.5 + 1);
%         
%         for bin_idx = bin_L : bin_H
%             LP_spec_L(:,bin_idx) = lpc(spec_L(bin_idx,:)',LP_spec_order); 
%             LP_spec_R(:,bin_idx) = lpc(spec_R(bin_idx,:)',LP_spec_order);
%         end
%         
%         
%         spec_recon_L = zeros(FFT_len*0.5+1, spec_recon_len);
%         spec_recon_R = zeros(FFT_len*0.5+1, spec_recon_len);
%         waveform_stack_L = zeros(FFT_len, spec_recon_len);
%         waveform_stack_R = zeros(FFT_len, spec_recon_len);
%         window_scale = FFT_len / (hop_size*2);
%         waveform_L = zeros(size_fshift+FFT_len,1);
%         waveform_R = zeros(size_fshift+FFT_len,1);        
%         
%         for j = 1 : spec_recon_len
%             for bin_idx = bin_L : bin_H
%                 for lp = 2:LP_spec_order
%                     %Left to right reconstruction
%                     spec_recon_L(bin_idx,j) = spec_recon_L(bin_idx,j) -1 * LP_spec_L(lp,bin_idx) * spec_L(bin_idx, spec_buff_len+j-(lp-1));
% 
%                     %Right to left reconstruction
%                     spec_recon_R(bin_idx,j) = spec_recon_R(bin_idx,j) -1 * LP_spec_R(lp,bin_idx) * spec_R(bin_idx, spec_buff_len+j-(lp-1));
%                 end
%             end
%             spec_recon_L(:,j) = abs(spec_recon_L(:,j));
%             spec_recon_R(:,j) = abs(spec_recon_R(:,j));
%             spec_L(:,spec_buff_len+j) = spec_recon_L(:,j);
%             spec_R(:,spec_buff_len+j) = spec_recon_R(:,j);
%             
%             waveform_stack_L(:,j) = real(ifft([spec_recon_L(:,j); flipud(spec_recon_L(2:FFT_len*0.5,j))]));
%             waveform_stack_R(:,j) = real(ifft([spec_recon_R(:,j); flipud(spec_recon_R(2:FFT_len*0.5,j))]));
%             
%             waveform_L(1+j*hop_size:FFT_len+j*hop_size) = waveform_L(1+j*hop_size:FFT_len+j*hop_size) +  waveform_stack_L(:,j);
%             waveform_R(1+j*hop_size:FFT_len+j*hop_size) = waveform_R(1+j*hop_size:FFT_len+j*hop_size) +  waveform_stack_R(:,j);
%         end
%         
%         waveform_L = waveform_L(FFT_len*0.5 + 1 : FFT_len*0.5 + size_fshift);
%         waveform_R = waveform_R(FFT_len*0.5 + 1 : FFT_len*0.5 + size_fshift);
%         s_proc_bp_rec = [waveform_L; waveform_R];%./ window_scale;
        

        s_proc_bs=train_bsf(s_frame_1d, fs, n, Wn);
        s_proc = s_proc_bp_rec + s_proc_bs;
        s_proc = s_proc .* w;
       
        %s_proc = zeros(size_frame,1);
        
        s_frame_rect_2d = s_frame_rect_1d;
    end
    
    %% remove tonal noise using notch filter
    if n_flag_prv == 1 || (notch_cnt > 1 && notch_cnt < TONE_len)
        s_proc = train_bsf(s_proc, fs, n, [5950 6150]); 
        notch_cnt = notch_cnt+1;
    end
    
    if notch_cnt == TONE_len
        notch_cnt = 0;
    end

    s_out(size_crnt : size_crnt-1 + size_frame) = s_out(size_crnt : size_crnt-1 + size_frame) + s_proc;
    
    
    
    %% Proccessed signal buffering
    s_frame_1d = s_frame;
    
    s_frame_rect_1d = s_frame_rect;
    dres_bp_1d = dres_bp;
    detect_mask_1d = detect_mask;
    n_flag_prv = n_flag;
    
    size_crnt = size_crnt + size_fshift;
    i = i+1;
end

wavwrite(s_out./32767, fs, 16, output);

% plot(s);
% hold on;
% plot(s_out(size_fshift:length(s_out))*100,'r');

end