global p;

%Codec settings
p.BC_GOLAY = 1;
p.enc_BPF = 1;
p.dec_BPF = 1;
p.dec_PowScale = 0;
p.INTERLEAVE = 1;
p.PLAY_PACKETWISE = 1;
p.MIC_ON = 1;
p.URL_Mode = 2; %1: Web, 2: Menu, Card
p.LOCALE = 'eng'; %Locale codes
% Tb = 120;
p.fs = 48000;
p.AMP = 15000;
p.init_frame = 100; %Initial frames for noise update
p.scale_thr = 1.2;
p.fc = 20000;
p.Rp = 1; %Packet rate: Databit / Packetbit: 2/3
p.Rc = 12/23; %Coding Rate: Mx12 in -> Mx23 out -> Golay coding : 12/23
% p.Rc = 1;
p.BPS_num = p.fs * p.Rp * p.Rc;

%Decoder setting
p.SCALE_MAX = 15000;
p.DEBUG = 1;

%%Parameters
p.fc_L = 18000;
p.fc_H = p.fs/2;
p.analysis_len = p.fs * 30; 
p.g = 2.0; %Spreading factor
p.B_init = (p.fc_L:p.fc_H ); %Initial target band, Hz
p.B_width_Unit = 5; %Unit half bandwidth to search optimum bandwidth, Bin
% p.A0 = 0.5; %Unit-normalizing constant for Attenuation factor
p.dist_scale = 1; %Scale distances
p.C_scale = 40964 / p.Rc;
p.N_guard_len = 5550;
p.N_sync_len = 0; %power preservation region that include sync (0 : off)
p.sil_len = 5000;

%% Attenuation model fitting
p.WORD_LEN = 50; %Maximum length of the unit word length
p.comp_len = 8192;
p.R_dist = 2; %Distance scaling by interger multiple
p.rand_in = 0;
p.test_iter = 100;
