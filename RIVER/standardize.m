function [X,m,st] = standardize(X)

% Standardize by column
% X: PxN matrix [P: variable(e.g. gene), N: instances (e.g. observations)] 

[n p] = size(X);

m = nanmean(X,2);
X = X - repmat(m,1,p);
% m = mean(X,2);
% X(isnan(X)) = 0;
% st = sqrt(sum(X.^2,2));
st = nanstd(X,0,2);
% X = X./sqrt(sum(X.^2,2)*ones(1,p));
X = X./repmat(st,1,p);