% This function returns local maxima of maxima lines.

function mxma_of_mxma_lines=find_maxima_of_maxima_lines(mxma_lin_cnnctd)
    
    ls=length(mxma_lin_cnnctd);
    mxma_of_mxma_lines=cell(1,3);
    
    for i=1:ls
        nb_lines=length(mxma_lin_cnnctd{i});
        mxma_of_mxma_lines{i}=struct;
        
        mxma_of_mxma_lines{i}.x=[];
        mxma_of_mxma_lines{i}.y=[];
        mxma_of_mxma_lines{i}.modulus=[];
        mxma_of_mxma_lines{i}.arguments=[];
%         mxma_of_mxma_lines{i}.scales=[];
        
        for j=1:nb_lines
            line_length=length(mxma_lin_cnnctd{i}{j}.modulus);
            [peaks,locs]=findpeaks(movmean(mxma_lin_cnnctd{i}{j}.modulus,5),...
                'MinPeakDistance',ceil(line_length/2));
            
            mxma_of_mxma_lines{i}.x=[mxma_of_mxma_lines{i}.x,mxma_lin_cnnctd{i}{j}.x(locs)'];
            mxma_of_mxma_lines{i}.y=[mxma_of_mxma_lines{i}.y,mxma_lin_cnnctd{i}{j}.y(locs)'];
            mxma_of_mxma_lines{i}.modulus=[mxma_of_mxma_lines{i}.modulus,peaks'];
            mxma_of_mxma_lines{i}.arguments=[mxma_of_mxma_lines{i}.arguments,...
                mxma_lin_cnnctd{i}{j}.arguments(locs)'];
%             mxma_of_mxma_lines{i}.scales=[mxma_of_mxma_lines{i}.scales,mxma_lin_cnnctd{i}{j}.z(locs)];
            
        end
    end
    
end