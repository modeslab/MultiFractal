% This function chains maxima of maxima lines through different scales. 

function chains=chain_mxma_of_mxma_lines_v2(mxma_of_mxma_lines,s)
    ls=length(mxma_of_mxma_lines);
    
    
    
    nchains=length(mxma_of_mxma_lines{1}.x);
    chains=cell(1,nchains);
    
    for j=1:nchains
        
        chains{j}=struct;
        
        chains{j}.x=mxma_of_mxma_lines{1}.x(j);
        chains{j}.y=mxma_of_mxma_lines{1}.y(j);
        chains{j}.arguments=mxma_of_mxma_lines{1}.arguments(j);
        chains{j}.modulus=mxma_of_mxma_lines{1}.modulus(j);
        chains{j}.scale=s(1);
        
    end
    
    survived_chains_idx=1:nchains;
    
    for i=1:ls-1
        disp(['scale ',num2str(i)])
        max_dist=s(i);
        [neighs_idx,dists]=knnsearch([mxma_of_mxma_lines{i}.x(:),mxma_of_mxma_lines{i}.y(:)],[mxma_of_mxma_lines{i+1}.x(:),mxma_of_mxma_lines{i+1}.y(:)]);
        
%         figure
%         scatter(mxma_of_mxma_lines{i}.x(:),mxma_of_mxma_lines{i}.y(:),[],1:length(mxma_of_mxma_lines{i}.x(:)),'filled')
%         hold on
%         scatter(mxma_of_mxma_lines{i+1}.x(:),mxma_of_mxma_lines{i+1}.y(:),[],neighs_idx,'d','filled')
%         
        rep_idx=[];
        
        for j=1:length(neighs_idx)
            tmp_idx=find(neighs_idx==neighs_idx(j));
            
            for k=1:length(tmp_idx)
                if (dists(j)~=min(dists(tmp_idx)))
                    rep_idx=[rep_idx,j];
                end
                if (dists(j)>max_dist)
                    rep_idx=[rep_idx,j];
                end

            end

        end
        
        neighs_idx(rep_idx)=[];
        mxma_of_mxma_lines{i+1}.y(rep_idx)=[];
        mxma_of_mxma_lines{i+1}.x(rep_idx)=[];
        mxma_of_mxma_lines{i+1}.modulus(rep_idx)=[];
        mxma_of_mxma_lines{i+1}.arguments(rep_idx)=[];
        
        survived_chains_idx=survived_chains_idx(neighs_idx);
        
        for j=1:length(neighs_idx)
%             plot([chains{survived_chains_idx(j)}.x(end),mxma_of_mxma_lines{i+1}.x(j)],[chains{survived_chains_idx(j)}.y(end),mxma_of_mxma_lines{i+1}.y(j)],'b-')
            
            chains{survived_chains_idx(j)}.x=[chains{survived_chains_idx(j)}.x,mxma_of_mxma_lines{i+1}.x(j)];
            chains{survived_chains_idx(j)}.y=[chains{survived_chains_idx(j)}.y,mxma_of_mxma_lines{i+1}.y(j)];
            chains{survived_chains_idx(j)}.modulus=[chains{survived_chains_idx(j)}.modulus,mxma_of_mxma_lines{i+1}.modulus(j)];
            chains{survived_chains_idx(j)}.arguments=[chains{survived_chains_idx(j)}.arguments,mxma_of_mxma_lines{i+1}.arguments(j)];
            chains{survived_chains_idx(j)}.scale=[chains{survived_chains_idx(j)}.scale,s(i+1)];

        end
        
        
    end
    
    
end