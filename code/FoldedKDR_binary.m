clearvars;

% data
n=200;
pl=5;
pr=5;
dl=1;
dr=1;
MU0=zeros(pl*pr,1);
SIGMA0=diag(ones(pl*pr,1));
L_true=[1;-1;1;0;zeros(pl-4,1)];
R_true=[-1;1;-1;0;zeros(pr-4,1)];
B_true =kron(R_true,L_true);

% Random initialization of projection matrix 
rand('state',1);
seed=floor(rand(1)*100);
randn('state', seed);

M=10;                     % Number of Iterations
s=10;                     % Number of Iterations of Starting Values                                         
mintr=1e+006;


for i=1:M
    
for j=1:n
    vecX(j,:)=mvnrnd(MU0,SIGMA0);
    X=reshape(vecX(j,:),pl,pr);
    if 0.1*(L_true'*X*R_true)^2+0.5>1
    Y(j,:)=1;
    else Y(j,:)=0;
    end
end

vecX=(vecX-repmat(mean(vecX),n,1))./repmat(std(vecX),n,1);
sgx0=2*(dl*dr)/(pl*pr)*MedianDist(vecX)^2;
sgy0=2*MedianDist(Y)^2;

for j=1:s

[L,R,B,tr]=FoldedKDR(vecX,Y,pl,dl,pr,dr,sgx0,sgy0);

if tr < mintr
   mintr=tr;
   minL =L;
   minR =R;
   minB =B;
end

L =minL;
R =minR;
B =minB;


end

vc1(i)=vctc(L,L_true);
vc2(i)=vctc(R,R_true);
LZC(i)=FoldedKDRLZC(B,B_true);

end

vc_mean=[mean(vc1),mean(vc2)]
LZC_mean=mean(LZC)
