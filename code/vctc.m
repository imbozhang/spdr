
function [vc]=vctc(A,B)
[A,D,V]=svd(A,0);   % Normalization
[B,D,V]=svd(B,0);   % Normalization
mat=B'*A*A'*B;
d  =eig(mat);
d  =(d+abs(d))/2;
vc =sqrt(prod(d));
tc =sqrt(mean(d));