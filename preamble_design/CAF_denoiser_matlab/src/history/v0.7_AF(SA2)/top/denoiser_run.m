clear;
clc;

mkdir('input');
mkdir('output');
prm.input = 'input/samsung/input/CAF_mix_48m.wav';
prm.output = 'output/CAF_mix_48m_matlab.wav';
prm.FS = 48000;
prm.SIZE_FRAME = 1024;
prm.SIZE_FSHIFT = 1024;
prm.LPorder = 6;
prm.LPorder_recon = 6;
prm.BP_cutL = 2500;
prm.BP_cutGAP = 500;
prm.BP_cutH = 12000;  
prm.MEAN_thr = 9;
prm.TONE_len = 2;
prm.RECON_POWER = 1.0;

prm.DRES_OUT =0;
prm.DETECT = 0;
prm.BYPASS = 0;
prm.BYPASS_INV = 0;
prm.TONE_REMOVE = 1;

addpath('src');

DRES_denoiser(prm);
