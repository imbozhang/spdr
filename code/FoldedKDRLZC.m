% LZC distance 
% Li,Zha and Chiaromonte (2005) 
function LZC1 = FoldedKDRLZC (A,B)

PA=A*pinv(A'*A)*A';
PB=B*pinv(B'*B)*B';
LZC1=sum(sum((PA-PB).*(PA-PB),1),2);
%LZC2=norm(PA-PB,'fro');