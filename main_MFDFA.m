clear
close all
clc

q=linspace(-3,3,50); % vector q

%pth='the_path_including_images_to_analyze';
list=dir; % list of files in the directory
nf=numel(list)-2; % number of files in the directory 

n_freq=2; % No. of frequencies to be removed in order to remove oscillatory trends

for i=1:nf
    
    fln=list(i+2).name; % filename to be processed
    disp(['step: ',num2str(i),' / ',num2str(nf),' -> ',fln])
    
    A=double(importdata(fln)); % reading the image 
    A=(A-mean(mean(A)))/std2(A); % normalization of the data
    
    sz=size(A);
    
    s1=floor(sqrt(sz(1)*sz(2)/1600)); % the lower limit of s
    s2=floor(sqrt(sz(1)*sz(2)/37)); % the upper limit of s
    scale=linspace(log2(s1),log2(s2),100);
    s=round(2.^(scale));
    s=unique(s); % vector s
%     
    
    X1=kill_low_freqs(A,n_freq); % removing the low frequencies in the data to get rid of oscillatory trends
    [Hq1,Fq1,err1]=MFDFA_2d_function(X1,s,q); % performing MFDFA
    [al1,f_al1,al_err1,f_al_err1]=sing_spec(q,Hq1,[Hq1-err1(1,:);err1(2,:)-Hq1]); % caculating the singulatory spectrum and corresponding errors
    clear X1
    
    figure
    errorbar(q,Hq1,Hq1-err1(1,:),err1(2,:)-Hq1)

%     X2=shuffle_2d(kill_low_freqs(A,n_freq)); % shuffling the data for null hypothesis 
%     [Hq2,Fq2,err2]=MFDFA_2d_function(X2,s,q); % performing MFDFA on the shuffled data
%     [al2,f_al2,al_err2,f_al_err2]=sing_spec(q,Hq2,[Hq2-err2(1,:);err2(2,:)-Hq2]); % caculating the singulatory spectrum and corresponding errors
%     clear X2
    
    
    X3=Rnkw_gaus(kill_low_freqs(A,n_freq)); % creating the ranked-wise Gaussian surrogate
    [Hq3,Fq3,err3]=MFDFA_2d_function(X3,s,q); % performing MFDFA on the surrogate
    [al3,f_al3,al_err3,f_al_err3]=sing_spec(q,Hq3,[Hq3-err3(1,:);err3(2,:)-Hq3]);
    clear X3
    
    
%     X4=RD_phase_surrogate(kill_low_freqs(A,n_freq)); % creating the the random-phased surrogate of the data
%     [Hq4,Fq4,err4]=MFDFA_2d_function(X4,s,q); performing MFDFA on the surrogate
%     [al4,f_al4,al_err4,f_al_err4]=sing_spec(q,Hq4,[Hq4-err4(1,:);err4(2,:)-Hq4]);
%     clear X4
%     clear A
    fln=erase(fln,".jpg");
    save([fln,'_res.mat']) % saving the analysis results
    
    figure('DefaultAxesFontSize',18)
    errorbar(al1,f_al1,f_al_err1(1,:),f_al_err1(2,:),al_err1(1,:),al_err1(2,:),'bo','DisplayName','Original Data')
    hold on
%    errorbar(al2,f_al2,f_al_err2(1,:),f_al_err2(2,:),al_err2(1,:),al_err2(2,:),'ro')
    errorbar(al3,f_al3,f_al_err3(1,:),f_al_err3(2,:),al_err3(1,:),al_err3(2,:),'go','DisplayName','RW Gaussian')
%    errorbar(al4,f_al4,f_al_err4(1,:),f_al_err4(2,:),al_err4(1,:),al_err4(2,:),'mo')
    xlabel('\alpha')
    ylabel('f(\alpha)')
    legend('Location','north east')
    title(fln)
end

load handel
sound(y,Fs)
