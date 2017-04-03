global p;

%Codec settings
p.BC_GOLAY = 2; %0 = off, 1 = Matlab golay, 2 = C_DLL Golay 
p.enc_BPF = 1;
p.dec_BPF = 1;
p.dec_PowScale = 0;
p.INTERLEAVE = 0;
p.PLAY_PACKETWISE = 0;
p.MIC_ON = 0;
p.URL_Mode = 1; %1: Web, 2: Menu, Card
p.LOCALE = 'eng'; %Locale codes
p.fs = 44100;
p.AMP = 20000;
p.init_frame = 100; %Initial frames for noise update
p.scale_thr = 1.0;
p.fc = 19600;
p.Tb = 810; %multiple of 9 
p.Rp = 1; %Packet rate: Databit / Packetbit: 2/3
p.Rc = 12/23; %Coding Rate: Mx12 in -> Mx23 out -> Golay coding : 12/23
% p.Rc = 1;
p.BPS_num = p.fs * p.Rp * p.Rc;
p.chirp_sync = 1; %Off: sequence-based sync, On: signal-based sync

%Decoder setting
p.SCALE_MAX =25000;
p.DEBUG = 0;

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

p.DRES_blk = 9;
p.t1 = 0.5*p.DRES_blk*p.Tb / p.fs; %~166ms %Half length of the preamble 
p.fl = 18050; %Chirp bound (low)
p.fh = 22050; %Chirp bound (high)
p.rr = 0.9; %Rectification ratio
p.ws = 1.0; % > 1 thinner, <1 : thicker
p.PTHR_MAIN = 1;
p.PTHR_RAKE = 5;
p.BETA = 2; %Power factor for correlated result
p.MP_MIN = 5; %Minimum duration betwen paths (in sample)


%% Attenuation model fitting
p.WORD_LEN = 50; %Maximum length of the unit word length
p.comp_len = 8192;
p.R_dist = 2; %Distance scaling by interger multiple
p.rand_in = 0;
p.test_iter = 100;
