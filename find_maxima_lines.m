% this function computes the maxima points in the Modulus field along
% direction of angles. nn describes the number of neighbors that are taken
% into account for comparison. In a syemmetric 

function [skel_tmp]=find_maxima_lines(Modulus,Angles,nn)



    szx=size(Modulus);
    skel_tmp=zeros(szx(1),szx(2));
    
    parfor xx=1+nn:szx(1)-nn
        
        tmp_v=zeros(1,szx(2));
        
        if(mod((xx-nn),100)==1)
            disp([num2str(xx-nn),'/',num2str(szx(1)-2*nn)])
        end
        for yy=1+nn:szx(2)-nn
            
            
            theta=Angles(xx,yy);
            xlen=nn*cos(theta);
            ylen=nn*sin(theta);
            xin_tmp=round(linspace(xx-xlen,xx+xlen,2*nn+1));
            yin_tmp=round(linspace(yy-ylen,yy+ylen,2*nn+1));
            lind=xin_tmp+(yin_tmp-1)*szx(1);
            M_tmp=Modulus(lind);
            
            
            mx_tst=islocalmax(M_tmp);
            
            if(mx_tst(1+nn))
                tmp_v(1,yy)=Modulus(xx,yy);
                
            end
        end
        skel_tmp(xx,:)=tmp_v(1,:);
        
    end
    
    


end