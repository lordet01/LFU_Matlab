function [flag, delay]=env_detect(s_buff, s_preamble, Tb)

H = length(s_preamble) / Tb;
flag = 0;
delay = 0;
r_max = 10;
r_LR = 0.1;

[corr_val,lag]=xcorr(s_buff, s_preamble);
corr_val = abs(corr_val + sqrt(-1)*hilbert(corr_val));
[max_val,idx]=max(corr_val);
mean_val_L = mean(corr_val(1:(H-2)*Tb));
mean_val_R = mean(corr_val((H)*Tb : end));

if max_val > r_max * mean_val_L && max_val > r_max * mean_val_R && abs(mean_val_R-mean_val_L) < r_LR * (mean_val_R+mean_val_L) &&...
        (idx > (H-2)*Tb && idx < (H)*Tb)
    
    delay_corr = lag(idx) + Tb;
    if delay_corr <= 0
        flag = 1;
        
        %Centering windowed symbol
        e_frame = abs(s_buff(1:Tb) + sqrt(-1)*hilbert(s_buff(1:Tb)));
        [~,delay_frame] = min(e_frame);
        if delay_frame > 0.5*Tb
            delay_frame = Tb - delay_frame;
        end
        delay = delay_corr + delay_frame;
        if delay > 0
            delay = 0;
        end
    end
end

end