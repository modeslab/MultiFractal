function [al,f_al]=sing_spec_from_tauq(q,tauq)
    lq=length(tauq);
    dq=q(2)-q(1);
    t_dif=zeros(1,lq-4);
    for i=1:lq-4
        t_dif(i)=-1*(tauq(i+4)-8*tauq(i+3)+8*tauq(i+1)-tauq(i))/(12*dq);
    end
    q=q(3:end-2);
    tauq=tauq(3:end-2);
    al=t_dif;
    
    f_al=q.*t_dif-tauq;
end