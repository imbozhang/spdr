profile clear;
profile on;
clearvars;
t=cputime;

% set the seed
% s = RandStream('swb2712','Seed',606);
% RandStream.setDefaultStream(s);


% data
n=100;
pl=5;
pr=5;
dl=1;
dr=1;
MU0=zeros(pl*pr,1);
SIGMA0=diag(ones(pl*pr,1));
L_true=[1;-1;0;0;zeros(pl-4,1)];
R_true=[-1;1;0;0;zeros(pr-4,1)];
B_true =kron(R_true,L_true);

ANL=1;                   % Maximum value for anealing
ETA=10;                  % Range of the golden ratio search
M=1;                     % Number of Iterations
                     

for i=1:M
    
for j=1:n
    e=normrnd(0,1,1,1);
    vecX(j,:)=mvnrnd(MU0,SIGMA0);
    X=reshape(vecX(j,:),pl,pr);
    Y(j,:)=0.1*(L_true'*X*R_true+2)^2+0.5*e;
end

vecX=(vecX-repmat(mean(vecX),n,1))./repmat(std(vecX),n,1);

sgx0=2*MedianDist(vecX)^2;
sgy0=2*MedianDist(Y)^2;

L0=randn(pl,dl);          % L is initialized
R0=randn(pr,dr);          % R is initialized

[L,R,B,tr]=FoldedKDR(vecX,Y,pl,dl,pr,dr,sgx0,sgy0,ANL,ETA,L0,R0);

LZC(i)=FoldedKDRLZC(B,B_true);
vc1(i)=vctc(L,L_true);
vc2(i)=vctc(R,R_true);
end

e=cputime-t;
LZC_median=median(LZC)
vc_median=[median(vc1),median(vc2)]
plot(tr,'bp');

%profile viewer;