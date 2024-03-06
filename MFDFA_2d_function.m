% This function performs MFDFA on the 2d data (image) X with scales s
% and orders q. The outputs are the generalized Hurst exponent Hq,
% dtrended fluctuation function Fq and err is the interval bounds of Hq 
% with .95 confidence level.

function [Hq,Fq,err]=MFDFA_2d_function(X,s,q)

ls=length(s); % number of scales
lq=length(q); % number of q's

Fq=zeros(lq,ls); 
q2=repmat(q,[ls,1]); % repeated q for daster access when run in parallel
parfor ns=1:ls
%     clc
%     disp(['itteration: ',num2str(ns),' / ',num2str(ls)])
    q_temp=q2(ns,:);
    
    segments=floor(size(X)/s(ns)); 
    segments1=segments(1); % number of segments in the first dimension
    segments2=segments(2); % number of segments in the second dimension
    
    RMS=zeros(segments1,segments2); % Root mean square of detrended fluctuations for each segment 
    for v=1:segments1
        for u=1:segments2
            Index1=((((v-1)*s(ns))+1):(v*s(ns))); % indeces of datapoints in the current segment 
            Index2=((((u-1)*s(ns))+1):(u*s(ns)));
            
            tmp=cumsum(X(Index1,Index2),1);
            tmp=cumsum(tmp,2);
            
            [XOut, YOut, ZOut] = prepareSurfaceData(Index1, Index2, tmp);
            
            [~,~,otpt]=fit([XOut,YOut],ZOut,'poly11'); % ploynomial fit for detrending 
            RMS(v,u)=sqrt(mean((otpt.residuals).^2)); 
        end
    end
    for nq=1:lq
        Fq(nq,ns)=mean(RMS(:).^q_temp(nq)).^(1/q_temp(nq));
    end

end
Hq=zeros(1,lq);
err=zeros(2,lq);
scale1=s;

for nq=1:lq
    fitobj=fit(log(scale1)',log(Fq(nq,:))','poly1');
    Hq(nq)=fitobj.p1;
    ci = confint(fitobj,0.95);
    err(:,nq)=ci(:,1);
end


end