function [V,D,XX,Q] = fun_TRCA_Matrix(X)

nChans  = size(X,1);
nTrials = size(X,3);
S = zeros(nChans, nChans);


X1 = X(:,:);
XX=sum(X,3);
S2 = XX*XX'-X1*X1';
X1 = X1 - repmat(mean(X1,2),1,size(X1,2));
Q = X1*X1';
[V_raw,D_raw] = eig(S2,Q);
eigvalue=diag(D_raw);
[D,index]=sort(eigvalue(:,1),1,'descend');
V=V_raw(:,index);
W=orth(V);
end