% This funcion produces Random-Phased data-surrogate. This is done by
% randomizing the phses in the Fourier transform and then taking inverse
% transform.

function Y=RD_phase_surrogate(X)
    X=X-mean(X(:));
    szx=size(X);
%     Y=zeros(size(X));
    r=exp(2*pi*1i*rand(szx));
    Y=ifft2(fft2(X).*r,'symmetric');
end