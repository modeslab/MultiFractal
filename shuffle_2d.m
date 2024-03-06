% This functions suffles the image to remove all correlations and keep the
% distribution preserved.

function Y=shuffle_2d(X)

rind=randperm(numel(X));

Y=reshape(X(rind),size(X));

end