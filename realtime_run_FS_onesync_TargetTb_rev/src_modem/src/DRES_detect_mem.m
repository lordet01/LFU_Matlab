function [g, prm]=DRES_detect_mem(p)


prm.SIZE_FRAME = floor(p.t1*p.fs * 2);
prm.SIZE_FSHIFT = floor(p.t1*p.fs * 2);
prm.LPorder = 10;
prm.MEAN_thr = p.PTHR_MAIN;


%Temporal buffers for operation
g.s_frame_1d = zeros(prm.SIZE_FRAME,1);
g.s_frame_LP_1d = zeros(prm.SIZE_FRAME*2,1);
g.peak_val_1d = 0;
g.dres_bp_1d = zeros(prm.SIZE_FRAME*2,1);

g.frame_idx = 1;
g.cnt = 1;
g.w = hanning(prm.SIZE_FRAME*2);
g.peak_pos = 1;
g.peak_val_sel = 0;

end