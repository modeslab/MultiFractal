% This function removes the lowest m frequencies in the data to get rid of
% the oscillatory trends. It works by setting the amplitude of the Fourier 
% transform to zero.

function fltrd_im=kill_low_freqs(X,m)

szx=size(X);
ft=fft2(X);
    

ft(1:m,1:m)=0;
ft(szx(1)-m+1:szx(1),1:m)=0;

ft(1:m,szx(2)-m+1:szx(2))=0;
ft(szx(1)-m+1:szx(1),szx(2)-m+1:szx(2))=0;

fltrd_im=ifft2(ft,'symmetric');