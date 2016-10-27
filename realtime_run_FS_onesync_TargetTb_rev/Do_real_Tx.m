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

while (1) %repeat send
    %% Cleanup outputs
    delete('src_modem/snd/*.wav');
    delete('src_modem/preamble/*.mat');
    delete('src_modem/output/*.txt');
    delete('src_modem/output/ref/*.txt');
    delete('src_modem/bit/*.txt');
    delete('src_modem/bit/ref/*.txt');
    delete('src_modem/rcv_sim_fs1/*.wav');
    delete('src_modem/rcv_sim_fs1/Tb_Test/*.wav');
    
    %% Generate arbitrary input txt
    if p.rand_in
        fin_r = fopen('src_modem/input/in_rnd.txt', 'w');

        symbols = ['a':'z' 'A':'Z' '0':'9'];
        MAX_ST_LENGTH = p.WORD_LEN;
        %         stLength = randi(MAX_ST_LENGTH);
        stLength = MAX_ST_LENGTH;
        nums = randi(numel(symbols),[1 stLength]);
        st = symbols (nums);
        fprintf(fin_r, st);
        fclose(fin_r);
        fname = 'in_rnd';
    else
        fname = 'menu';
    end

    fc = p.fc;
    Tb = p.Tb;
    %% Generate Tx signal
    fname_txt = [fname,'.txt'];
    fname_snd = [fname];
    for j = 1:length(fc)
        for i = 1:length(Tb)
            if Tb(i) < 100
                fname_snd_full = [fname_snd,'_Tb',num2str(Tb(i)),'_Fc',num2str(fc(j))];
            else
                fname_snd_full = [fname_snd,'_Tb',num2str(Tb(i)),'_Fc',num2str(fc(j))];
            end
            fexist = fopen(['src_modem/snd/',fname_snd_full,'.wav']);
            if fexist == -1
                u_mod(fname_txt, [fname_snd_full,'.pcm'], Tb(i), fc(j), p);
            end
        end
    end
    fclose all;
    
    pause(1);
    if p.PLAY_PACKETWISE == 0
        break;
    end
end
