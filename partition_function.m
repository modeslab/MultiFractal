% This function calculates the partition function Z of a given skeleton
% named chains. It also fits a power-law to each curve corresponding to q
% in order to calculate tau_q. change the range of fit in the marked line
% to get more accurate result. It is recommended to use MFDFA for better
% estimate of multifractality. 


function [tau_q,err,Z]=partition_function(chains,s,q)

Z=zeros(length(s),length(q));
tau_q=zeros(size(q));
err=zeros(2,length(q));
for is=1:length(s)

    for k=1:length(chains)

        if ismember(s(is),chains{k}.scale)
            for iq=1:length(q)
                Z(is,iq)=Z(is,iq)+...
                    (max(chains{k}.modulus(chains{k}.scale<=s(is))))^q(iq);
            end
        end
    end

end

for iq=1:length(q)
    tmp_idx=find((s>2^5.2).*(s<2^5.7)); % Change the range of fit here.
    
    fitobj=fit(log(s(tmp_idx)),log(Z(tmp_idx,iq)),'poly1');
    % fitobj=fit(log(s((Z(:,iq)~=0))),log(Z((Z(:,iq)~=0),iq)),'poly1');
    tau_q(iq)=fitobj.p1;
    ci = confint(fitobj,0.99);
    err(1,iq)=tau_q(iq)-ci(1,1);
    err(2,iq)=tau_q(iq)-ci(2,1);
end





end