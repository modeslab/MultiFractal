% This peice of code is the same as "observables.m" but it prints the
% results into a .txt file.

clear 
clc
close all

list2=dir;

nf2=(numel(list2)-2)/2;
% orig_d_alpha=zeros(1,nf2);
% orig_d_alpha_int=orig_d_alpha;
% RK_d_alpha=zeros(1,nf2);
% RK_d_alpha_int=RK_d_alpha;

fileID = fopen('/Users/bahadori/Documents/sci_project/MFDFA/new_sorted_data/BAPN/BAPN_res.txt','a');

% fprintf(fileID,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n','age','lit_no','whole width',... % uncomment this line bfore runing this code for the first time for each stage to prduce titles of the columns
%     'left side range','whole integral','left side integral','max loc',...
%     'max alpha','min alpha','file name');

age=15.5; % change this 
lit_no=1;

for cntr=1:nf2
    flnm=list2(cntr*2+2).name;
    load(flnm)
    
    
    [tmp,mx_idx]=findpeaks(f_al1);
%     figure
%     plot(al1,f_al1,'.-',al1(mx_idx),tmp,'ro')
%     hold on
    if(length(mx_idx)>1)
        
        dists=min(mx_idx-1,length(al1)-mx_idx);
        
        mx_idx=mx_idx(dists==max(dists));
    end
%     plot(al1(mx_idx),2,'g^')
    
    orig_d_alpha=max(al1);

    orig_d_alpha_whole_int=abs(sum(f_al1(1:end-1).*(diff(al1(1:end)))));
    orig_d_alpha_left_int=abs(sum(f_al1(mx_idx:end-1).*(diff(al1(mx_idx:end)))));

    
    fprintf(fileID,'%.2f\t %.0f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %s\n',age, lit_no, max(al1)-min(al1),al1(mx_idx)-min(al1),orig_d_alpha_whole_int,orig_d_alpha_left_int,al1(mx_idx),max(al1),min(al1),flnm);
    
    fprintf('%.2f\t %.0f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %s\n',age, lit_no, max(al1)-min(al1),al1(mx_idx)-min(al1),orig_d_alpha_whole_int,orig_d_alpha_left_int,al1(mx_idx),max(al1),min(al1),flnm);
    
     
end





fclose(fileID);



% figure(1)
% scatter(age*ones(1,nf2),orig_d_alpha,'bo')
% hold on
% scatter(age*ones(1,nf2),RK_d_alpha,'ro')
% 
% figure(2)
% scatter(age*ones(1,nf2),orig_d_alpha_int,'bo')
% hold on
% scatter(age*ones(1,nf2),RK_d_alpha_int,'ro')