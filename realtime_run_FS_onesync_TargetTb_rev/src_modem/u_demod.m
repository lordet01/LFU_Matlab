function u_demod(fname_rcv, fname_txt, bit_name_txt, Tb, p)

addpath('src_modem/src');
load('src/BPF_44100_19200-20200.mat');

%load golay code libary from C
addpath('src_modem/lib');
if not(libisloaded('golay_lib'))
    loadlibrary('golay_lib.dll','golay.h')
end
%libfunctions('golay_lib')

%calllib('golay_lib', 'decode_golay', 53);


MIC_ON = p.MIC_ON;
INTERLEAVE = p.INTERLEAVE;
BC_GOLAY = p.BC_GOLAY;
Fs = p.fs;
SCALE_MAX = p.SCALE_MAX;

%Parameters for Golay coding
M = 2; %2 x 12 bit
N = 3; %3 x 8 bit
B = 8; %Unit symbol length
C = 12; %Unit code length
G = 23; %Unit golay code length

    
if BC_GOLAY
    G = G;
else
    G = C;
end

header_bits = [1 1 1 0 0 1 1 0];
ender_bits = [1 1 1 1 1 1 1 0];
bit_buff = zeros(1, M*G);
bit_blk_c = zeros(M, C);

s_prv = zeros(1, Tb);
n_avg = 0;

s_BPF_tail_prv_1 = zeros(1, length(BPF)-1);
header_flag = 0;
err_cnt_char = 0;
bit_cnt = 0;

fin = fopen(fname_rcv, 'rb');
fout = fopen(fname_txt, 'wt');
fbit = fopen(bit_name_txt, 'wt');
fprintf(fbit, '\n');

if p.DEBUG
    fdebug1 = fopen('src_modem/snd/rcv.pcm', 'ab');
    fdebug2 = fopen('src_modem/snd/rcv_Proc.pcm', 'ab');
    fdebug3 = fopen('src_modem/snd/rcv_Proc_Scale.pcm', 'ab');
end


if MIC_ON
    %Initialize Audio recording
    [obj]=dsp_record(0, Fs, 1, Tb, 'start');
end


%Write buffer (to check ender for the text
c_buff = zeros(p.WORD_LEN,1); %Maximum length for one sentence to be sent
c_idx = 0;

% Discard wav header 
if MIC_ON == 0
    fread(fin, 22, 'int16');
end

nFrmIdx = 0;
url_cnt = 1;
url_buff = zeros(p.WORD_LEN,1);
while (1)
    if MIC_ON
        [obj, s]=dsp_record(obj, Fs, 1, Tb, 'record');
        s = s .* 32767;
        s = s';
        len = length(s);
    else
        %File-based processing
        [s, len] = fread(fin, Tb, 'int16');
        s = s';
    end
    if p.DEBUG
        fwrite (fdebug1, s, 'int16');
    end    
    %Check Ultrasonic band's average power
    if len ~= Tb
        break;
    else
        if p.dec_BPF
            s_BPF = conv(s,BPF);
            s_BPF(1:length(s_BPF_tail_prv_1)) = s_BPF(1:length(s_BPF_tail_prv_1)) + s_BPF_tail_prv_1;
            s_BPF_tail_1 = s_BPF(length(s)+1:length(s_BPF));
            
            s_BPF = s_BPF(1:length(s));
            
            if p.dec_PowScale
                %Noise estimation by averaging frames
                if nFrmIdx < p.init_frame
                    n_avg = n_avg + sum(s_BPF.^2);        
                elseif nFrmIdx == p.init_frame
                    n_avg = n_avg / p.init_frame;
                else
                    %frame-wise noise tracking
    %                 if SNR_ind == 0
    %                     n_avg = eta * n_avg + (1-eta) * s_BPF.^2;
    %                 end 
                end

                %Calculate A-posteriori SNR
    %             disp('sig');
    %             disp(sum(s_BPF.^2));
    %             disp('noise');
    %             disp(sum(n_avg));
                if sum(s_BPF.^2) > p.scale_thr*n_avg
                    SNR_ind = 1;
                else
                    SNR_ind = 0;
                end

                %Write preprocessed waveform
                if p.DEBUG
                    fwrite (fdebug2, s_BPF, 'int16');
                end
                if (SCALE_MAX > 0) && SNR_ind
                    s = SCALE_MAX .* (s_BPF / max(abs(s_BPF)));
                end
% %                 s = s_BPF;
% %                 s_BPF = conv(s,BPF);
% %                 s_BPF(1:length(s_BPF_tail_prv_2)) = s_BPF(1:length(s_BPF_tail_prv_2)) + s_BPF_tail_prv_2;
% %                 s_BPF_tail_2 = s_BPF(length(s)+1:length(s_BPF));
% % 
% %                 s_BPF = s_BPF(1:length(s));
            end
            s = s_BPF;
            
            %Write preprocessed waveform
            if p.DEBUG
                fwrite (fdebug3, s, 'int16');
            end
        end
    end
    
% %     S = fft(s, Tb);
% %     bin_B = p.fs / Tb;
% %     B_idx = floor(20000 / bin_B);
% %     tmp = S(B_idx-1:B_idx+1);
% %     S = zeros(1,Tb/2);
% %     S(B_idx-1:B_idx+1) = tmp;
% %     S_proc = [0, S, flipud(S(1:end-1)')'];
% %     s = real(ifft(S_proc, Tb));

%     disp([s(6), s(11), s_prv(6), s_prv(11)]);
    bit = bit_detect(s_prv, s);

    %Stack bits
    bit_buff(1:M*G-1) = bit_buff(2:M*G);
    bit_buff(M*G) = bit;
    bit_cnt = bit_cnt + 1;
    
    if header_flag == 0 || (header_flag == 1 && bit_cnt == M*G)
        %Pefrom channel decoding
        bit_blk_ibc = matB2C(bit_buff, M, G);
        
        if INTERLEAVE
            bit_blk_ibc = [bit_blk_ibc'; [0, 0] ]'; %Zero padding
            bit_blk_bc = matdeintrlv(bit_blk_ibc',12,2);
            bit_blk_bc = bit_blk_bc';
            bit_blk_bc = bit_blk_bc(:,1:23);
        else
            bit_blk_bc = bit_blk_ibc;
        end
        
        
        if BC_GOLAY == 1
            %[23, 12] Golay decoding
            bit_blk_c(1,:) = golaycodec(bit_blk_bc(1,:));
            bit_blk_c(2,:) = golaycodec(bit_blk_bc(2,:));
        elseif BC_GOLAY == 2
           deci_tmp1 = bi2de(bit_blk_bc(1,:), 'left-msb'); 
           deci_tmp2 = bi2de(bit_blk_bc(2,:), 'left-msb');
           deci_out1 = calllib('golay_lib', 'decode_golay', deci_tmp1);     
           deci_out2 = calllib('golay_lib', 'decode_golay', deci_tmp2);     
           
           bit_blk_c(1,:) = de2bi(deci_out1,C,'left-msb');
           bit_blk_c(2,:) = de2bi(deci_out2,C,'left-msb');
        else
            bit_blk_c = bit_blk_bc;
        end
        bit_blk = matB2C(bit_blk_c, N, B);        
        bit_cnt = 0;

        if header_flag
            bit_write(bit_blk, fbit);
            for k = 1 : N
                if mod(sum(bit_blk(k,:)),2)==0;
                    c = char(bi2de(bit_blk(k,1:7),'left-msb'));
                    fprintf('%c', c);
                    c_idx = c_idx + 1;
                    c_buff(c_idx) = c;
                    url_buff(url_cnt) = c;
                    url_cnt = url_cnt + 1;
                else
                    err_cnt_char = err_cnt_char + 1;
                end
            end
            %Finish decoding (Successful transmission)
            if isequal(ender_bits, bit_blk(3,:)) && c_idx <= p.WORD_LEN
                if p.URL_Mode == 1 %Web
                    URL = char(url_buff);
                    [~,URL]=strtok(URL,'@');
                    [URL,~]=strtok(URL,10);
                    flag = strncmp(URL,'http',4);
                    if flag == 0
                        try
                            URL = ['https://',URL'];
                        catch
                            URL = ['https://',URL];
                        end
                    end
                    fprintf('\n%s\n',URL);
%                     web(URL, '-new', '-notoolbar');
                    system(['start chrome ','"',URL,'"']);
                    
                    url_buff = zeros(p.WORD_LEN,1);
                    url_cnt = 1;
                elseif p.URL_Mode == 2 %Menu or Card
                    URL = char(url_buff);
                    [~,URL]=strtok(URL,'@');
                    [URL,~]=strtok(URL,10);
                    URL = [URL',p.LOCALE,'.pdf'];
                    flag = strncmp(URL,'http',4);
                    if flag == 0
                        try
                            URL = ['https://',URL'];
                        catch
                            URL = ['https://',URL];
                        end
                    end
                    disp(['\n',URL]);
                    
                    [key1,~]=size(strfind(URL, 'http'));
                    [key2,~]=size(strfind(URL, 'pdf'));
                    if key1 == 1 && key2 == 1
                        system(['start chrome ','"',URL,'"']);
                    end
                    url_buff = zeros(p.WORD_LEN,1);
                    url_cnt = 1;
                end
                
                header_flag = 0;
                c_idx = 1;
                fprintf(fout, '\n');
                
                fprintf(fbit, '\n');
            end
            
            %Finish decoding (Force termination)
            if c_idx > p.WORD_LEN %Forcely terminate current
                fprintf(fout, 'FAILED');
                header_flag = 0;
                c_idx = 1;
                fprintf(fout, '\n');
            end
        else
            %Sync set
            if isequal(header_bits, bit_blk(1,:)) && isequal(header_bits, bit_blk(2,:))&& isequal(header_bits, bit_blk(3,:))
                header_flag = 1;
                    
                %Kill & start player
                [~, ~]= system('taskkill /F /IM chrome.exe');
            end
        end
    end
%     disp([bi2de(bit_blk(1,:),'left-msb'),bi2de(bit_blk(2,:),'left-msb'),bi2de(bit_blk(3,:),'left-msb')]);
    
    
    s_prv = s;
    if p.dec_BPF
        s_BPF_tail_prv_1 = s_BPF_tail_1;

    end
    
    if nFrmIdx <= p.init_frame
        nFrmIdx = nFrmIdx + 1;
    end
end
fclose('all');
unloadlibrary('golay_lib');

end
