function y=train_bpf(x, Fs, n, Wn)

Fn = Fs/2;

[b, a] = butter(n, Wn/Fn);

y=filter(b,a,x);


end