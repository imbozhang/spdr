clearvars;

% Random initialization of projection matrix 
rand('state',1);
seed=floor(rand(1)*100);
randn('state', seed);

Data=dlmread('EEG.data');

pl=256;
pr=64;
dl=2;
dr=1;
vecX=Data(:,2:(pl*pr+1));      
Y=Data(:,1);
n=length(Y);

%% pca or not
pca=1;                         

if pca==1
    
    pl=30;
    pr=20;

load('pcavecX.mat')
    
end


% Normalization

vecX=(vecX-repmat(mean(vecX),n,1))./repmat(std(vecX),n,1);


%% 5 fold cross validation or not

Nrep=100;

for i=1:Nrep
    
   index=randsample(n,n);
   
   for j=1:5
       index_test=index(25*(j-1)+1:min(25*j,n));
       index_train=index;
       index_train(index_test)=[];
       n_test=length(index_test);
       n_train=length(index_train);
       
       vecX_train=vecX(index_train,:);
       vecX_test=vecX(index_test,:);
       Y_train=Y(index_train);
       Y_test=Y(index_test);
       
       
       % key part
       
       sgx0=2*(dl*dr)/(pl*pr)*MedianDist(vecX_train)^2;
       sgy0=2*MedianDist(Y_train)^2;
       
       [L,R,B,tr]=FoldedKDR(vecX_train,Y_train,pl,dl,pr,dr,sgx0,sgy0);
       
       Z_train=vecX_train*B;
       
       %RGN=max([max(abs(Z_train(:,1))), max(abs(Z_train(:,2)))])*1.2;
       %fxy = figure('Position',[0 400 300 300]);
       %axis([-RGN RGN -RGN RGN]);
       %axis square;
       %KDRplot(fxy,vecX_train,Y_train,B);
       
       X11_train = Z_train(:,1);
       X12_train = Z_train(:,2);
       group = Y_train;
       Z_test=vecX_test*B;
       X11_test = Z_test(:,1);
       X12_test = Z_test(:,2);       
       
       C = classify([X11_test X12_test],[X11_train X12_train],group,'quadratic');
       
       err(i,j)=length(find(C-Y_test))/n_test;
       
       [i,j,err(i,j)]
 
   end
   
end
       
       err_kdr=mean(mean(err))
  


