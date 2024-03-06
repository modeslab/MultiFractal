clear
clc
close all

%% producing the image from brownian function
H=.6;
disp('generating brownian field...')
[field1,~,~,~]=Brownian_field(H,2^13);

A=field1(1:2500,1:2500);

clear field1 
%% producing the image with one singularity and one local maximum

% res=1025*2-1;
% A=zeros(res);
% clc
% x=-floor(res/2):floor(res/2); % linspace(-500,500,res);
% y=x;
% for i=1:res
%     for j=1:res
%         A(i,j)=exp(-((x(i)-256)^2+(y(j)-256)^2)/(2*128^2))-sqrt((x(i)+256)^2+(y(j)+256)^2)^.3;
%     end
% end


%% loading data from filename

% filename='2021_10_28_E17,5_B008_#1_right.jpg';
% A=double(importdata(filename));
% A=randn(3e3,3e3);

%% initializing

X=A;
szx=size(X)-[500,500];%

clear A

scale=linspace(log2(30),log2(250),100);
s=(2.^(scale));
s=unique(s(:));


%% transform calculation
disp('wavelet transform calculation...')
cwtstruct = cwtft2(X,'wavelet',{'mexh',{2,1,1}},'scales',s,'norm','L2');

% an is the argument field
an=zeros(szx(1),szx(2),length(s));
% M is the modulus field
M=an;

fx=an;
fy=an;


% to use the whole image
% [fx(:,:,:),fy(:,:,:)]=gradient(cwtstruct.cfs(:,:,1,:));
% an(:,:,:)=angle(-1j*fx+fy);
% M(:,:,:)=abs(fx+1j*fy);


% to cut the margins of the image
[fx(:,:,:),fy(:,:,:)]=gradient(cwtstruct.cfs(251:end-250,251:end-250,1,:));
an(:,:,:)=angle(1j*fx+fy);
M(:,:,:)=abs(fx+1j*fy);

figure
imagesc(X(251:end-250,251:end-250))
axis xy;
colormap gray

% plot an example of each field.
figure
subplot(2,2,1);
imagesc(fx(:,:,1)); axis xy;
colormap gray
subplot(2,2,2);
imagesc(fy(:,:,1)); axis xy;
colormap gray
subplot(2,2,3);
imagesc(M(:,:,1)); axis xy;
colormap gray
hold on
subplot(2,2,4);
imagesc((an(:,:,1))); axis xy;
colormap gray

clear cwtstruct

WTMM=zeros(szx(1),szx(2),length(s)); % Wavelet Transform Modulus Maxima






% nn is number of neighbours from each side that the point is maximum among them 
nn=2;



%% 
disp('Maxima line finding...')

Modulus(:,:)=zeros(szx(1),szx(2));
Args(:,:)=zeros(szx(1),szx(2));

for is=1:length(s)
    tic
    disp(['step: ',num2str(is),'/',num2str(length(s))])
    
    Modulus(:,:)=M(:,:,is);
    Args(:,:)=an(:,:,is);
    
    skel_tmp=find_maxima_lines(Modulus,Args,nn);

    WTMM(:,:,is) =skel_tmp(:,:);
    toc
end




%%
disp('making structure of maxima lines')
[maxima_line_struct]=make_struct_of_maxima_lines(WTMM,an,s);

clear M an
save('tst_mexh')

%%
disp('connecting maxima lines...')
tic
mxma_lin_cnnctd=connect_maxima_lines_v3(maxima_line_struct);
toc

%%
disp('finding maxima of the connected maxima lines...')
mxma_of_mxma_lines=find_maxima_of_maxima_lines(mxma_lin_cnnctd);


%% ploting the maxima lines and their maxima
disp('ploting the maxima lines and their maxima...')


minZ=log(min(s));
figure
%imagesc(X);
%colormap gray
%imagesc(M(:,:,1));
hold on
%plot3(x_ind,y_ind,log(z_ind)-min(log(z_ind)),'ro','MarkerSize',1)
colormap copper
for is=1:length(s)
    number_of_lines=length(mxma_lin_cnnctd{is});
    for iq=1:number_of_lines
        tmp_inds=find((mxma_lin_cnnctd{is}{iq}.x<1500).*(mxma_lin_cnnctd{is}{iq}.y<1500));
        plot3(mxma_lin_cnnctd{is}{iq}.x(tmp_inds),mxma_lin_cnnctd{is}{iq}.y(tmp_inds),log(s(is)*ones(size(tmp_inds)))-minZ)%,5,log(mxma_lin_cnnctd{i}{j}.modulus),'filled')
    end
    tmp_inds=find((mxma_of_mxma_lines{is}.x<1500).*(mxma_of_mxma_lines{is}.y<1500));    
    scatter3(mxma_of_mxma_lines{is}.x(tmp_inds),mxma_of_mxma_lines{is}.y(tmp_inds),log(s(is)*ones(size(tmp_inds)))-minZ,15,'r','filled')
    
    
end
colorbar



%% chaining the maximas of maxima lines
disp('chaining the maximas of maxima lines and plotting...')

chains=chain_mxma_of_mxma_lines_v2(mxma_of_mxma_lines(1,1:50),s(1:50));


disp('plotting skeleton...')
xrng=[0,1500];
yrng=[0,1500];

minZ=log(min(s));
figure
imagesc(X); axis xy;
colormap gray
% imagesc(M(:,:,1)); axis xy;
hold on
% colormap copper
for is=1:length(s)
    number_of_lines=length(mxma_lin_cnnctd{is});
%     for j=1:number_of_lines
%         tmp_inds=find((mxma_lin_cnnctd{i}{j}.x<xrng(2)).*(mxma_lin_cnnctd{i}{j}.y<yrng(2))...
%             .*(mxma_lin_cnnctd{i}{j}.x>xrng(1)).*(mxma_lin_cnnctd{i}{j}.y>yrng(1)));
%         plot3(mxma_lin_cnnctd{i}{j}.x(tmp_inds),mxma_lin_cnnctd{i}{j}.y(tmp_inds),log(s(i)*ones(size(tmp_inds)))-minZ)%,5,log(mxma_lin_cnnctd{i}{j}.modulus),'filled')
%     end
    tmp_inds=find((mxma_of_mxma_lines{is}.x>xrng(1)).*(mxma_of_mxma_lines{is}.y>yrng(1))...
        .*(mxma_of_mxma_lines{is}.x<xrng(2)).*(mxma_of_mxma_lines{is}.y<yrng(2)));    
    scatter3(mxma_of_mxma_lines{is}.x(tmp_inds),mxma_of_mxma_lines{is}.y(tmp_inds),log(s(is)*ones(size(tmp_inds)))-minZ,[],'r','filled');%,mxma_of_mxma_lines{is}.modulus(tmp_inds),)
    
    
end
% figure
for is=1:length(chains)
    if( length(chains{is}.scale)>=2 && chains{is}.x(1)<xrng(2) && chains{is}.y(1)<yrng(2) ...
            && chains{is}.x(1)>xrng(1) && chains{is}.y(1)>yrng(1))
    plot3(chains{is}.x,chains{is}.y,log(chains{is}.scale)-minZ,'bo-')
    hold on
    end
    
end
ylim(yrng)
xlim(xrng)


disp('plotting power laws...')
x_sings=[];
y_sings=[];
exp_sings=[];
exp_err=[];
% figure('DefaultAxesFontSize',18)

for is=1:length(chains)
    
    if (length(chains{is}.scale)>=4)
        p=fit(log(chains{is}.scale(1:4))',log(chains{is}.modulus(1:4))','poly1');
%         disp(num2str(i))
        ci = confint(p,0.95);
%         plot(log(chains{is}.scale),log(chains{is}.modulus),'.-','DisplayName',...
%             ['exponent: ',num2str(p.p1),'  ( ',num2str(ci(1,1)),'--',num2str(ci(2,1)),' )'])
%         hold on

        if (abs(p.p1)>3*abs(ci(1,1)-ci(2,1)))
            x_sings=[x_sings,chains{is}.x(1)];
            y_sings=[y_sings,chains{is}.y(1)];
            exp_sings=[exp_sings,p.p1];
            exp_err=[exp_err,[ci(1,1);ci(2,1)]];
        end
    end
    
end
% legend('Location','northwest')

disp('plotting the image along with i singularities...')

thrsld=-10;
figure
ax1 = axes;
imagesc(X(251:end-250,251:end-250)); axis xy;
colormap(ax1,'gray')
hold on
ax2 = axes;
scatter(x_sings(exp_sings>thrsld),y_sings(exp_sings>thrsld),[],...
    exp_sings(exp_sings>thrsld),'filled')
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);


%% Calculating partition function and multifractal measures
disp('Calculating partition function and multifractal measures...')
q=linspace(-5,5,40);

[tau_q,err_tauq,Z]=partition_function(chains,s,q);

Hq=(tau_q+2)./q;
err_hq=abs(err_tauq./(q));

[al,f_al]=sing_spec_from_tauq(q,tau_q);

clist = colormap(jet(length(q)));
figure('DefaultAxesFontSize',18)
for i=1:length(q)
plot(log2(s),log2(Z(:,i)),'.-','Color',clist(i,:))
hold on
end
colormap(clist)
c=colorbar('Ticks',[0,1],'TickLabels',[q(1),q(end)]);
c.Label.String = 'q';
xlabel('Log2(s)')
ylabel('Log2(Z(s,q))')





figure('DefaultAxesFontSize',18)
% tauq_2=Hq_2.*q -2;
% err_tauq2=abs(err_Hq_2.*q);
% errorbar(q,tauq_2,err_tauq2(1,:),err_tauq2(2,:),'bo','DisplayName','MFDFA')
hold on
errorbar(q,tau_q+2,err_tauq(1,:),err_tauq(2,:),'ro','DisplayName','WTMMMM')
plot(q,q*H-2,'g','DisplayName','Theory')
legend('Location','north west')
xlabel('q')
ylabel('\tau (q)')

figure('DefaultAxesFontSize',18)
errorbar(q,Hq,err_hq(1,:),err_hq(2,:))



figure('DefaultAxesFontSize',18)
plot(al,f_al,'bo-')
xlabel('\alpha')
ylabel('f(\alpha)')
% xlim([0,1])
% ylim([0,2])

%% compare with MFDFA

disp('comparing with MFDFA...')
s_MFDFA= unique(round(2.^(linspace(log2(64),log2(400),120))));
[Hq_2,Fq,err_Hq_2]=MFDFA_2d_function(X,s_MFDFA,q);
Hq_2=Hq_2-2;

figure('DefaultAxesFontSize',18)
loglog((s_MFDFA),(Fq))
xlabel('Log2(s)')
ylabel('Log2(F(q,s))')


figure('DefaultAxesFontSize',18)
errorbar(q,Hq_2,err_Hq_2(1,:),err_Hq_2(2,:))
xlabel('q')
ylabel('h(q)')



figure('DefaultAxesFontSize',18)
tauq_2=Hq_2.*q -2;
err_tauq2=abs(err_Hq_2.*q);
errorbar(q,tauq_2,err_tauq2(1,:),err_tauq2(2,:),'bo','DisplayName','MFDFA')
hold on
errorbar(q,tau_q+2,err_tauq(1,:),err_tauq(2,:),'ro','DisplayName','WTMMMM')
plot(q,q*H-2,'g','DisplayName','Theory')
legend('Location','north west')
xlabel('q')
ylabel('\tau (q)')

[al2,f_al2]=sing_spec_from_tauq(q,tauq_2);
figure('DefaultAxesFontSize',18)
plot(al2,f_al2,'bo-')
xlabel('\alpha')
ylabel('f(\alpha)')




%% saving final data
disp('saving...')
save('tst_mex')

load handel
sound(y,Fs)
