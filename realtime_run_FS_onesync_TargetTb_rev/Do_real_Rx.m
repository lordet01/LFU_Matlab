% /*********************************************
%      GIST LF-Ultrasonic Communication Channel Simulator (lvl1_FreeSpace Model) (v1.0)
%      Human Media Communication & Processing Laboratory (HuCom)
%   Developer: Kwang Myung Jeon
%     
%      Change Log.
%                     - [20150809](v1.0) Initial released
%                     - [20150809](v1.0) Goal: Get Eb/N0 for different
%                     distances(l) and SNRs
%                     - [20150809](v1.0) Refered article [1]
%                        [1] M. Stojanovic, "On the relationship between capacity and distance in
%                        an underwater acoustic communication channel," 2006.
%                     - [20150831](v1.1) Merged Ultrasonic modem 
%                     - [20150914](v1.2) Added Sync header on whole sequence, instead of each packet 
%                                 [H;D;D]...[H;D;D]...[H;D;D]  --->   [H;H;H]...[D;D;D]...[D;..E]
%                     - [20150915](v1.3) Revised Modem's synchronization conditions for improved 
%                                        robustness under various delay and background noise 
%                     - [20150923](v1.4) Added envelope detector and PN sequence based synchronization mode 
%                                 
% ***********************************************/
fclose('all');
close('all');
clear;
clc;
addpath('src');
addpath('src_modem');
initial_setting;


%% Decode RX signals ( Apply Tx signal on channel model)

snd_path = 'src_modem/snd';
rcv_path = 'src_modem/rcv_sim_fs1';
snd_list = dir(snd_path);
snd_list = snd_list(3:end);


Tb_target = p.Tb;
for i = 1:length(Tb_target)
    %Decode Rx signal
%     fname_wav_rcv = 'src_modem/rcv_sim_fs1/URL_out.wav';
    fname_wav_rcv = 'src_modem/snd/music_testTV_218_sample.wav';
    [~,fname_rcv] = strtok(fname_wav_rcv, '/');
    [~,fname_rcv] = strtok(fname_rcv, '/');
    [fname_rcv,~] = strtok(fname_rcv, '/');
%     fname_snd_name = strtok(fname_snd,'.');
    fname_txt_rcv = ['src_modem/output/',fname_rcv,'.txt'];
    bit_name_txt_rcv = ['src_modem/bit/',fname_rcv,'.txt'];
%     preamble_rcv = ['src_modem/preamble/p_',fname_snd_name,'.pcm.mat'];
    u_demod(fname_wav_rcv, fname_txt_rcv, bit_name_txt_rcv, Tb_target(i), p);
end


