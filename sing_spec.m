% This function computes the singularity spectrum and corresponding errors
% from generalized Hurst exponent. 

function [al,f_al,al_err,f_al_err]=sing_spec(q,Hq,Hq_err)
    lq=length(Hq);
    dq=q(2)-q(1);
    h_dif=zeros(1,lq-4);
    h_diff_err=zeros(2,lq-4);
    for i=1:lq-4
        h_dif(i)=-1*(Hq(i+4)-8*Hq(i+3)+8*Hq(i+1)-Hq(i))/(12*dq);
        h_diff_err(:,i)=-1*(Hq_err(:,i+4)-8*Hq_err(:,i+3)+8*Hq_err(:,i+1)-Hq_err(:,i))/(12*dq) + dq^4;
    end
    q=q(3:end-2);
    Hq=Hq(3:end-2);
    Hq_err=Hq_err(:,3:end-2);
    
    al= Hq+q.*h_dif;
    al_err=Hq_err+repmat(q,2,1).*h_diff_err;
    f_al= 2 + q.^2 .* h_dif;
    f_al_err=h_diff_err.*repmat(q,2,1).^2;
end