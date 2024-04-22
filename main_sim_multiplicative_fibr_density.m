% this code simulates ossification process for a fixed set of parameters and
% visualizes the process over iteration steps. 

clear
clc
close all

%% initializing

n_iter= 2000; % number of iterations for ossification process
im_size= [900,1500]; % size of the image.



beta_in=3e-4; % the inverse of length scale of phosphate gradientn beta
phos_grad= repmat(exp(- beta_in * (1:im_size(2))),[im_size(1),1]); % Exponential phospate gradient
figure
surf(phos_grad,'EdgeColor','none')
view(2)
set(gca,'DataAspectRatio',[1 1 1e-3])

colorbar


% col_fib=make_triangular_lattice(im_size(1),im_size(2),100);
%col_fib=rndm_ln(im_size(1),im_size(2),100,1.5e3); previous sims were with
%this

% alpha=1000;
% col_fib=(rndm_ln_ordr_grad(im_size(2),im_size(1),alpha,100,cll_fbr_num))';





nighborhood_cntr=4e-3; % the relative contribution of the local positive feedback.

crr_len=20; % the kernel width for the local positive feedback. 

fibr_neigh=12; % the kernel width for the effect of collagen fibers. 
epsilon=0; % a small offset for ossification probablity to add white-noise if needed


col_fib_cntr=3e-2; % the relative contribution of the collagen fiber field.
omega=1; % the power for the local feedback. Throughout our study this is set to 1 but can be altered to study nonlinear feedbacks


s1=floor(sqrt(im_size(1)*im_size(2)/2000));
s2=floor(sqrt(im_size(1)*im_size(2)/200));
scale=linspace(log2(s1),log2(s2),100);
s=round(2.^(scale));
s=unique(s); % scales for MFDFA


q=linspace(-3,3,50); % The orders for MFDFA
%% simulation

n_rel=4; % number of realizations 

cll_fbr_num=1000:200:3400; % number of collagen fibers
cll_fbr_len=16:2:34; % length of collagen fibers


cl_no_no=length(cll_fbr_num);
cl_len_no=length(cll_fbr_len);

col_fib_all=cell(cl_no_no,cl_len_no,n_rel);
im_all=cell(cl_no_no,cl_len_no,n_rel);


Hq=cell(cl_no_no,cl_len_no,n_rel);
Fq=cell(cl_no_no,cl_len_no,n_rel);
err=cell(cl_no_no,cl_len_no,n_rel);

for i=1:cl_no_no
    for j=1:cl_len_no
        
        parfor k=1:n_rel
            
            disp(['i=',num2str(i),'/',num2str(cl_no_no),...
                '   j=',num2str(j),'/',num2str(cl_len_no), ...
                '   k=',num2str(k),'/',num2str(n_rel)])
            
            col_fib=rndm_ln(im_size(1),im_size(2),cll_fbr_len(j),cll_fbr_num(i));
            col_fib_all{i,j,k}=col_fib;
%             figure
%             imagesc(col_fib) 
%             set(gca,'DataAspectRatio',[1 1 1])
            im= zeros(im_size);
            P_fibr=col_fib_cntr*smoothdata(smoothdata(col_fib,1,'movmean',fibr_neigh),2,'movmean',fibr_neigh);

            for t=1:n_iter

                prob_mat=(nighborhood_cntr*(smoothdata(smoothdata(im.^omega,1,'movmean',crr_len),2,'movmean',crr_len)) ...
                     +P_fibr)...
                    .* ( phos_grad);





                prob_mat=prob_mat/max(max(prob_mat))+epsilon;


                rnd_mat=rand(im_size);

                im=im + (rnd_mat<prob_mat);



            end
            
            figure
            imagesc(im)
            set(gca,'DataAspectRatio',[1 1 1])
            xticks([])
            yticks([])
            im_all{i,j,k}=im;
            
            
            try
                disp(['doing MFDFA','i=',num2str(i),'/',num2str(cl_no_no),...
                '   j=',num2str(j),'/',num2str(cl_len_no), ...
                '   k=',num2str(k),'/',num2str(n_rel)])
            
                [Hq{i,j,k},Fq{i,j,k},err{i,j,k}]=MFDFA_2d_function(im,s,q);
            catch
            end
            
            
        end
        save('vs_density_8rel')
    end
end












%% MFDFA measures
al=cell(cl_no_no,cl_len_no,n_rel);
f_al=cell(cl_no_no,cl_len_no,n_rel);
al_err=cell(cl_no_no,cl_len_no,n_rel);
f_al_err=cell(cl_no_no,cl_len_no,n_rel);

MF_mes=zeros(cl_no_no,cl_len_no,n_rel);
for i=1:cl_no_no
    for j=1:cl_len_no
        
        for k=1:n_rel    
            try
            [al{i,j,k},f_al{i,j,k},al_err{i,j,k},f_al_err{i,j,k}]=...
                sing_spec(q,Hq{i,j,k},[Hq{i,j,k}-err{i,j,k}(1,:);...
                err{i,j,k}(2,:)-Hq{i,j,k}]);

            
            MF_mes(i,j,k)=max(al{i,j,k})-min(al{i,j,k});
            catch
            end
        end
    end
end



figure
heatmap(cll_fbr_num,cll_fbr_len,mean(MF_mes,3))

% figure
% errorbar(q,Hq,Hq-err(1,:),err(2,:)-Hq)
% 
% figure('DefaultAxesFontSize',18)
% errorbar(al,f_al,f_al_err(1,:),f_al_err(2,:),al_err(1,:),al_err(2,:),'bo','DisplayName','Original Data')
% 
% 
% figure
% loglog(s,Fq)

%% saving
save('vs_density_8rel')