function data=sub_linear(data)
%ignore constant and linear component
A=[ones(size(data,2),1),(1:size(data,2))'];
[A,s]=svd(A,0);
ind=eye(size(data,2));
%ignore constant and linear component
data=data-data*A*A';
% data=data*(ind-A*A');
