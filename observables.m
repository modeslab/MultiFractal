% this code can be run in the directory of the data after main_MFDFA.m is
% done runing. This will calculate the measures such as spectrum width for
% all *_res.mat files.

clear 
clc

list2=dir;

nf2=(numel(list2)-3)/2;
orig_d_alpha=zeros(1,nf2);
orig_d_alpha_int=orig_d_alpha;
RK_d_alpha=zeros(1,nf2);
RK_d_alpha_int=RK_d_alpha;
for cntr=1:nf2
    flnm=list2(cntr*2+3).name;
    load(flnm)
    
    mx_idx=find(f_al1==max(f_al1));
    
%     orig_d_alpha(cntr)=al1(mx_idx)-min(al1);
%     RK_d_alpha(cntr)=al3(mx_idx)-min(al3);
%     orig_d_alpha(cntr)=al1(mx_idx);
%     RK_d_alpha(cntr)=al3(mx_idx);
    orig_d_alpha(cntr)=max(al1);
    RK_d_alpha(cntr)=max(al3);

    
%     orig_d_alpha_int(cntr)=sum(f_al1(mx_idx:end-1).*(diff(al1(mx_idx:end))));
%     RK_d_alpha_int(cntr)=sum(f_al3(mx_idx:end-1).*(diff(al3(mx_idx:end))));

    orig_d_alpha_int(cntr)=min(al1);
    RK_d_alpha_int(cntr)=min(al3);

end

age=17.5;

figure(1)
scatter(age*ones(1,nf2),orig_d_alpha,'bo')
hold on
scatter(age*ones(1,nf2),RK_d_alpha,'ro')

figure(2)
scatter(age*ones(1,nf2),orig_d_alpha_int,'bo')
hold on
scatter(age*ones(1,nf2),RK_d_alpha_int,'ro')