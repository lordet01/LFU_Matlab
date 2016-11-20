function y=train_bsf(x, Fs, n, Wn)

Fn = Fs/2;

[b, a] = butter(n, Wn/Fn, 'stop');

y=filter(b,a,x);


end