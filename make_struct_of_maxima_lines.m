% This function makes a structure out of WTMM for the further analyses.

function [maxima_line_struct]=make_struct_of_maxima_lines(WTMM,args,s)

    szx=size(WTMM);
    
    if (length(szx)<3)
        szx(3)=1;
    end
    
    maxima_line_struct=cell(1,szx(3));

    
    for i=1:szx(3)

        maxima_line_struct{i}=struct;


        inds=find(WTMM(:,:,i)~=0);
        [I,J] = ind2sub(szx,inds);
        m_tmp=WTMM(:,:,i);
        arg_tmp= args(:,:,i);


        maxima_line_struct{i}.y_inds=I;
        maxima_line_struct{i}.x_inds=J;
        maxima_line_struct{i}.modulus=m_tmp(inds);
        maxima_line_struct{i}.arguments=arg_tmp(inds);
        maxima_line_struct{i}.z=(s(i)*ones(length(inds),1));



    end

end