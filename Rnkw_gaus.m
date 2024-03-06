% This function produces the Ranked-wise Gaussian data-surrogate. The data
% is sorted first and each value is then replaced by that of a sorted
% random array with Gaussian distribution.

function Y=Rnkw_gaus(X)

[~,ind]=sort(X(:));

r=randn(numel(X),1)*std(X(:));
rsrt=sort(r(:));
Y_tmp=zeros(numel(X),1);
Y_tmp(ind(:)) = rsrt;
Y=reshape(Y_tmp, size(X,1), size(X,2));
end