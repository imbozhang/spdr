clearvars;

% data
n=200;
pl=10;
pr=10;
dl=2;
dr=2;
mu=1;
sigma2=0.1;
tao2=1.5;
MU0=zeros(pl*pr,1);
SIGMA0=diag([1;sigma2;ones(pr-2,1);sigma2;ones((pl-1)*pr-1,1)]);
MU1=[mu;zeros(pl,1);mu;zeros(pl*(pr-1)-2,1)];
SIGMA1=diag([1;tao2;ones(pr-2,1);tao2;ones((pl-1)*pr-1,1)]);
e1=[1;0;zeros(pl-2,1)];
e2=[0;1;zeros(pr-2,1)];
L_true=[[1;0;0;0;zeros(pl-4,1)],[0;1;0;0;zeros(pl-4,1)]];
R_true=[[1;0;0;0;zeros(pr-4,1)],[0;1;0;0;zeros(pr-4,1)]];
B_true =[kron(e1,e1),kron(e1,e2),kron(e2,e1),kron(e2,e2)];

% Random initialization of projection matrix 
rand('state',1);
seed=floor(rand(1)*100);
randn('state', seed);

M=1;                      % Number of Iterations
s=3;                      % Number of Iterations of Starting Values                                         


for i=1:M
    
mintr=1e+006;

Y = binornd(1,0.5,n,1);

for j=1:n
    if Y(j)==0
       vecX(j,:) = mvnrnd(MU0,SIGMA0);
    else
       vecX(j,:) = mvnrnd(MU1,SIGMA1);
    end
end



vecX=(vecX-repmat(mean(vecX),n,1))./repmat(std(vecX),n,1);

sgx0=2*(dl*dr)/(pl*pr)*MedianDist(vecX)^2;
sgy0=2*MedianDist(Y)^2;

for j=1:s

[L,R,B,tr]=FoldedKDR(vecX,Y,pl,dl,pr,dr,sgx0,sgy0);

if tr< mintr
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

vc_LZC_mean=[mean(vc1),mean(vc2),mean(LZC)]