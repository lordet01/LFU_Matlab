function [y]=fconv(x, h)
%FCONV Fast Convolution
%   [y] = FCONV(x, h) convolves x and h.  The output of this 
%         function is scaled.
%
%      x = input vector
%      h = input vector
% 
%      See also CONV
%
%Version 2.0
%2003-2004, Stephen G. McGovern

Ly=length(x)+length(h)-1;   % 
Ly2=pow2(nextpow2(Ly));     % Find smallest power of 2 that is > Ly
% m=max(abs(x));      
X=fft(x, Ly2);		    % Fast Fourier transform
H=fft(h, Ly2);		    % Fast Fourier transform
Y=X.*H;        	           
y=real(ifft(Y, Ly2));       % Inverse fast Fourier transform
y=y(1:1:Ly);                % Take just the first N elements
% m=m/max(abs(y));
% y=m*y;
