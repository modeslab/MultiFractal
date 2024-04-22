% this code simulates ossification process for a fixed set of parameters and
% visualizes the process over iteration steps. 

clear
clc
close all

%% initializing

n_iter= 2000; % number of iterations for ossification process
im_size= [900,1500]; % size of the image.

im= zeros(im_size);


beta_in=3e-4; % the inverse of length scale of phosphate gradientn beta
phos_grad= repmat(exp(- beta_in * (1:im_size(2))),[im_size(1),1]); % Exponential phospate gradient
figure
surf(phos_grad,'EdgeColor','none')
view(2)
set(gca,'DataAspectRatio',[1 1 1e-3])

colorbar

cll_fbr_num=600; % number of collagen fibers.
col_fib=rndm_ln(im_size(1),im_size(2),100/5,2.5*cll_fbr_num); % randomly distributing collagen fibers with uniform distribution.

% col_fib=make_triangular_lattice(im_size(1),im_size(2),100); % creating a triangular mesh for collagen fiber field


% gamma=1000; % the length scale of the gradient for collagen fiber density
% col_fib=(rndm_ln_ordr_grad(im_size(2),im_size(1),gamma,100,cll_fbr_num))'; % Creat an exponential distribution for collagen fiber field 


figure
imagesc(col_fib) 
set(gca,'DataAspectRatio',[1 1 1])
col_fib_cntr=3e-2;



nighborhood_cntr=4e-3; % the relative contribution of the local positive feedback.

crr_len=20; % the kernel width for the local positive feedback. 
fibr_neigh=12; % the kernel width for the effect of collagen fibers. 
epsilon=0; % a small offset for ossification probablity to add white-noise if needed

omega=1; % the power for the local feedback. Throughout our study this is set to 1 but can be altered to study nonlinear feedbacks
%% simulation

% creat a video for the ossification simulation.
% vidfile = VideoWriter('40gaus_kernel_lin_corr_only_grad_mutiplied.mp4','MPEG-4');
% open(vidfile);
% figure

P_fibr=col_fib_cntr*smoothdata(smoothdata(col_fib,1,'movmean',fibr_neigh),2,'movmean',fibr_neigh); % The probability field corresponding to collagen fibers

for t=1:n_iter
    
    prob_mat=(nighborhood_cntr*(smoothdata(smoothdata(im.^omega,1,'movmean',crr_len),2,'movmean',crr_len)) ...
         +P_fibr)...
        .* ( phos_grad);



        
    
    prob_mat=prob_mat/max(max(prob_mat))+epsilon;
    
    
    rnd_mat=rand(im_size);
    
    im=im + (rnd_mat<prob_mat);
    
    
%     imagesc(im)
%     set(gca,'DataAspectRatio',[1 1 1])
%     xticks([])
%     yticks([])
%     title(['iteration:',num2str(t)])
%     
%     %     colormap('gray')
%     colorbar
%     caxis([0 1200])
%     frame = getframe(gcf);
% 
%     writeVideo(vidfile, frame);

    if (mod(t,100)==0)
        clc
        disp(['t= ',num2str(t),'/',num2str(n_iter)])
        
    end
end

% close(vidfile)


figure
imagesc(im)
set(gca,'DataAspectRatio',[1 1 1])
xticks([])
yticks([])
% im=abs(im+.1*mean2(im)*randn(size(im)));



%% MFDFA

    
s1=floor(sqrt(im_size(1)*im_size(2)/2000));
s2=floor(sqrt(im_size(1)*im_size(2)/200));
scale=linspace(log2(s1),log2(s2),100);
s=round(2.^(scale));
s=unique(s);


q=linspace(-3,3,50);

[Hq,Fq,err]=MFDFA_2d_function(im,s,q);
[al,f_al,al_err,f_al_err]=sing_spec(q,Hq,[Hq-err(1,:);err(2,:)-Hq]);

figure
errorbar(q,Hq,Hq-err(1,:),err(2,:)-Hq)

figure('DefaultAxesFontSize',18)
errorbar(al,f_al,f_al_err(1,:),f_al_err(2,:),al_err(1,:),al_err(2,:),'bo','DisplayName','Original Data')


figure
loglog(s,Fq)

%% saving
% save('alpha1000_corrlen20_collfiber600_omega1_3')