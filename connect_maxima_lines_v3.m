% This function connects WTMM's so that the local Maxima could be found
% along them. Uncomment the plotting lines to visualize the intermediate
% steps.

function mxma_lin_cnnctd=connect_maxima_lines_v3(maxima_line_struct)
    
    ls=length(maxima_line_struct);
    mxma_lin_cnnctd=cell(1,ls);
    
    min_len=6;
    
    
    
    parfor i =1:ls
        
        max_rad=.9*maxima_line_struct{i}.z(1);
        
%         figure
%         scatter3(maxima_line_struct{i}.x_inds,maxima_line_struct{i}.y_inds,...
%             maxima_line_struct{i}.arguments,[])
%         view(2)
        
       
        
        [neigh_idx,~]=rangesearch( [maxima_line_struct{i}.x_inds, maxima_line_struct{i}.y_inds]...
            ,[maxima_line_struct{i}.x_inds, maxima_line_struct{i}.y_inds],max_rad);
        
        
        init_len=length(maxima_line_struct{i}.x_inds);

        init_idx=1:init_len;
        

        while (sum(init_idx~=0)>2)
            clc
            disp(['scale: ',num2str(i),'/',num2str(ls)])
            disp([num2str(sum(init_idx~=0)),' / ',num2str(init_len)])
            
            remaining_idx=find(init_idx~=0);
            
            [orders]=...
                connect_points_v4( maxima_line_struct{i}.arguments,... 
                max_rad,remaining_idx(1),neigh_idx,maxima_line_struct{i}.x_inds,maxima_line_struct{i}.y_inds);

            init_idx(orders)=0; 
            
            if (length(orders)>min_len)

                

            

                mxma_lin_cnnctd{i}=[mxma_lin_cnnctd{i},cell(1,1)];

                mxma_lin_cnnctd{i}{end}=struct;

                x_cnctd=maxima_line_struct{i}.x_inds(orders);
                mxma_lin_cnnctd{i}{end}.x= x_cnctd;
                
                y_cnctd=maxima_line_struct{i}.y_inds(orders);
                mxma_lin_cnnctd{i}{end}.y=y_cnctd;
                
                M_cnctd=maxima_line_struct{i}.modulus(orders);
                mxma_lin_cnnctd{i}{end}.modulus=M_cnctd;
                
                an_cnctd=maxima_line_struct{i}.arguments(orders);
                mxma_lin_cnnctd{i}{end}.arguments=an_cnctd;
                
                mxma_lin_cnnctd{i}{end}.scales=...
                    maxima_line_struct{i}.z(1:length(x_cnctd));

                
                
                
                
            else
%                 plot(x_cnctd,y_cnctd,'go-')
                
                

            end
            
            
            
%             figure
%             plot(maxima_line_struct{i}.x_inds,maxima_line_struct{i}.y_inds,'bo')
%             hold on
%             plot(mxma_lin_cnnctd{i}{end}.x,mxma_lin_cnnctd{i}{end}.y,...
%                 'r.-')
            

            maxima_line_struct{i}.x_inds(orders)=-50;
            maxima_line_struct{i}.y_inds(orders)=-50;
            

        end

        
        
    end
    
end