function u_mod(fname_txt, fname_snd, Tb, fc, p)

addpath('src_modem/src');
% load('src/HPF.mat');
load('src/BPF_19500_20500.mat');

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
bit_blk = zeros(N,B);
blk_cnt = 0;
sync_on = 0;  %Header flag
sync_off = 0;
y_BPF_tail_prv = zeros(1, length(BPF)-1);

b_prev = 1;
header_bits = [1 1 1 0 0 1 1 0];
ender_bits = [1 1 1 1 1 1 1 0];
frame_idx = 0;
w = sqrt(hann(Tb, 'symmetric'));
% w = ones(Tb,1);
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


%% Ultrasonic generation. unit block = [header, char. char]
while (1)
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
            b_char = de2bi(d, 7);

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

        if BC_GOLAY
            bit_blk_bc = golaycodec(bit_blk_c);
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

end