function u_mod(fname_txt, fname_snd, Tb, fc, p)

addpath('src_modem/src');
% load('src/HPF.mat');
load('src/BPF_44100_19200-20200.mat');

%load golay code libary from C
if p.BC_GOLAY == 2
    addpath('src_modem/lib');
    if not(libisloaded('golay_lib'))
        loadlibrary('golay_lib.dll','golay.h')
    end
end
%libfunctions('golay_lib')

%calllib('golay_lib', 'decode_golay', 53);

BC_GOLAY = p.BC_GOLAY;
INTERLEAVE = p.INTERLEAVE;
PLAY_PACKETWISE = p.PLAY_PACKETWISE;
fs = p.fs;
AMP = p.AMP;
% Tb = 120;
% fc = 22200;

%Initialize Audio playing
[obj]=dsp_play(0, fs, Tb, 'start', 0);

%Parameters for Golay coding
M = 2; %2 x 12 bit 
N = 3; %3 x 8 bit
B = 8; %Unit symbol length
C = 12; %Unit code length
if BC_GOLAY > 0
    G = 23; %Golay bit length
else
    G = C;
end
bit_blk = zeros(N,B);
bit_blk_bc = zeros(M,G);
blk_cnt = 0;
sync_on = 0;  %Header flag
sync_off = 0;
y_BPF_tail_prv = zeros(1, length(BPF)-1);

b_prev = 1;
header_bits = [1 1 1 0 0 1 1 0];
ender_bits = [1 1 1 1 1 1 1 0];
frame_idx = 0;

%Window design
t_h = Tb;
ws = 1.0; % > 1 thinner, <1 : thicker
rr = 0.9; %Rectification ratio
t_win = floor(t_h*(1-rr));
w_tail = hann(t_win);
w_tail = w_tail .^ ws;
w = ones(t_h,1);
w(1:floor(t_win/2)) = w_tail(1:floor(t_win/2));
w(end-floor(t_win/2)+1:end) = w_tail(1+floor(t_win/2):2*floor(t_win/2));

% w = sqrt(hann(Tb, 'symmetric')); %Smoothing window
% w = ones(Tb,1); %Bypass windows
w = w';

% fname_txt = 'input_long.txt';
% fname_snd = 'snd_input_Tb120_Fc22.pcm';
fin = fopen(['src_modem/input/',fname_txt], 'rt');
fout = fopen(['src_modem/snd/',fname_snd], 'wb');
fbit = fopen(['src_modem/bit/bit_snd.txt'], 'wt');
fbit_bc = fopen(['src_modem/bit/bit_bc_snd.txt'], 'wt');

% %% Dummy codes (24bits)
% dummy_bits = [0 1 0 1 0 1 0 1 ...
%               0 1 0 1 0 1 0 1 ...
%               0 1 0 1 0 1 0 1];
% bit_blk_c = matB2C(dummy_bits, M, C);
% if BC_GOLAY
%     bit_blk_bc = golaycodec(bit_blk_c);
% else
%     bit_blk_bc = bit_blk_c;
% end
% [b_prev frame_idx]=encode_DPSK(bit_blk_bc, b_prev, frame_idx, fbit_bc, fout, w, Tb, fs, fc, AMP);

%Generate Chirp Sync preamble
if p.chirp_sync == 1
    t = -p.t1:1/fs:p.t1;
    y1 = chirp(t,p.fl,p.t1,p.fh,'quadratic',[],'concave');
    y2 = chirp(t,p.fl,p.t1,p.fh,'quadratic',[],'concave');

    t_h = floor(length(t) / 2);
    t_win = floor(t_h*(1-p.rr));
    w_tail = hann(t_win);
    w_tail = w_tail .^ p.ws;
    w_chirp = ones(t_h,1);
    w_chirp(1:floor(t_win/2)) = w_tail(1:floor(t_win/2));
    w_chirp(end-floor(t_win/2)+1:end) = w_tail(1+floor(t_win/2):2*floor(t_win/2));
    y1 = y1(1:t_h) .* w_chirp';
    y2 = y2(t_h+2:end) .* w_chirp';

    y_sync = [y1 , y2]';  
end


%% Ultrasonic generation. unit block = [header, char. char]
while (1)
    if p.chirp_sync == 1
        if sum(y_sync) ~= 0;
            fwrite(fout, floor(y_sync(1:p.Tb) .* AMP), 'int16');
            y_sync = [y_sync(p.Tb+1:end); zeros(p.Tb,1)];
            continue;
        end
    end
    
    if sync_on < N
        %Add header
        bit_blk(1:2, :) = bit_blk(2:3, :); %blk shift
        bit_blk(3, :) = header_bits;
        sync_on = sync_on + 1;
    else
        [content len] = fread(fin, 1);
        if len ~= 1
            %Add ender
            bit_blk(1:2, :) = bit_blk(2:3, :); %blk shift
            bit_blk(3, :) = ender_bits;
            sync_off = sync_off + 1;
        end
        
        %Add data
        if sync_off == 0
            d = double(content);
            b_char = de2bi(d, 7,'left-msb');

            cnt = 0;
            for k = 1:7
                if b_char(k) == 1
                    cnt = cnt + 1;
                end
            end
            if mod(cnt,2) == 1
                parity = 1;
            else
                parity = 0;
            end
            b_seq = [b_char parity];

            %Write data
            bit_blk(1:2, :) = bit_blk(2:3, :); %blk shift
            bit_blk(3, :) = b_seq; %blk update
        end
    end
    blk_cnt = blk_cnt + 1;
    
    %Perform Modulation
    if blk_cnt == N
        bit_write(bit_blk, fbit);
        bit_blk_c = matB2C(bit_blk, M, C);

        if BC_GOLAY == 1
            bit_blk_bc(1,:) = golaycodec(bit_blk_c(1,:));
            bit_blk_bc(2,:) = golaycodec(bit_blk_c(2,:));
        elseif BC_GOLAY == 2
           deci_tmp1 = bi2de(bit_blk_c(1,:), 'left-msb'); 
           deci_tmp2 = bi2de(bit_blk_c(2,:), 'left-msb');
           deci_out1 = calllib('golay_lib', 'encode_golay', deci_tmp1);     
           deci_out2 = calllib('golay_lib', 'encode_golay', deci_tmp2);     
           
           bit_blk_bc(1,:) = de2bi(deci_out1,G,'left-msb');
           bit_blk_bc(2,:) = de2bi(deci_out2,G,'left-msb');
        else
            bit_blk_bc = bit_blk_c;
        end

        if INTERLEAVE
            bit_blk_ibc = [bit_blk_bc'; [0, 0] ]'; %Zero padding
            bit_blk_ibc = matintrlv(bit_blk_ibc',12,2);
            bit_blk_ibc = bit_blk_ibc';
            bit_blk_ibc = bit_blk_ibc(:,1:23);
        else
            bit_blk_ibc = bit_blk_bc;
        end
        
        [y, b_prev frame_idx]=encode_DPSK(bit_blk_ibc, b_prev, frame_idx, fbit_bc, w, Tb, fs, fc, AMP);
        
        if p.enc_BPF
            y_BPF = conv(y,BPF);
            y_BPF(1:length(y_BPF_tail_prv)) = y_BPF(1:length(y_BPF_tail_prv)) + y_BPF_tail_prv;
            y_BPF_tail = y_BPF(length(y)+1:length(y_BPF));
            y_BPF = y_BPF(1:length(y));
            y = y_BPF;
            y_BPF_tail_prv = y_BPF_tail;
        end
        
        fwrite(fout, y, 'int16');
        if PLAY_PACKETWISE
            %y = y ./ 32767;
%             U = audioplayer(y, fs, 16);
%             playblocking(U);
            [obj]=dsp_play(obj, fs, Tb, 'record', y);
            disp('--send packet--');
        end
        
        blk_cnt = 0;
        
        if sync_off > N
            break;
        end
    end
end

%% Ultrasonic Sending (Filewise)
if PLAY_PACKETWISE == 0
    fout = fopen(['src_modem/snd/',fname_snd], 'rb');
    [Y] = fread(fout, inf, 'int16');
    %Y = Y ./ 32767;
    [obj]=dsp_play(obj, fs, Tb, 'record', Y');
    disp('--send packet--');
end

disp('---------Send Finished!---------');

%Convert PCM to WAV
pcm2wav(['src_modem/snd/',fname_snd],p);

% bit_blk(1:2, :) = bit_blk(2:3, :); %blk shift
% bit_blk(3, :) = ender_bits; %blk update
% blk_cnt = blk_cnt + 1;
% 
% while blk_cnt < N
%     bit_blk(1:2, :) = bit_blk(2:3, :); %blk shift
%     bit_blk(3, :) = zeros(1,B); %blk update
%     blk_cnt = blk_cnt + 1;
%     
%     %Perform Modulation
%     if blk_cnt == N
%         bit_write(bit_blk, fbit);
%         bit_blk_c = matB2C(bit_blk, M, C);
%         if BC_GOLAY
%             bit_blk_bc = golaycodec(bit_blk_c);
%         else
%             bit_blk_bc = bit_blk_c;
%         end
%         [b_prev frame_idx]=encode_DPSK(bit_blk_bc, b_prev, frame_idx, fbit_bc, fout, w, Tb, fs, fc, AMP);
%         blk_cnt = 0;
%         break;
%     end
% end

% 
% %% Mix signals for test purpose
% fEnvIn = fopen('mix/classic.pcm', 'rb');
% fEnvOut = fopen('snd/snd_url2_classic.pcm', 'wb');
% fout = fopen(['snd/',fname_snd], 'rb');
% delay = 10000; %Sample
% powlevel = 0.0; %power ratio
% 
% [Env] = fread(fEnvIn, inf, 'int16');
% [Stream] = fread(fout, inf, 'int16');
% 
% if length(Stream) > 0.5*length(Env)
%     while length(Env) < length(Stream)
%         Env = [Env; Env];
%     end
% end
% Stream = [zeros(delay,1); Stream];
% len_res = length(Env) - length(Stream);
% Stream = [Stream; zeros(len_res,1)];
% Mix = powlevel*(Env + Stream);
% fwrite(fEnvOut, Mix, 'int16');
% Mix = [Mix, zeros(size(Mix))];
% wavwrite(Mix./32767, 48000, 'snd/snd_classic.wav')

fclose('all');
unloadlibrary('golay_lib');

end