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
pca=0;                         

if pca==1
    
    pl=30;
    pr=20;

load('pcavecX.mat')
    
end

%% key part

vecX=(vecX-repmat(mean(vecX),n,1))./repmat(std(vecX),n,1);

sgx0=2*(dl*dr)/(pl*pr)*MedianDist(vecX)^2;
sgy0=2*MedianDist(Y)^2;

[L,R,B,tr]=FoldedKDR(vecX,Y,pl,dl,pr,dr,sgx0,sgy0);

Z=vecX*B;

RGN=max([max(abs(Z(:,1))), max(abs(Z(:,2)))])*1.2;
fxy = figure('Position',[0 400 300 300]);
axis([-RGN RGN -RGN RGN]);
axis square;
KDRplot(fxy,vecX,Y,B);

%% quadratic discriminant analysis

qda=1;

if qda==1

X11 = Z(:,1);
X12 = Z(:,2);
group = Y;
figure();
h1 = gscatter(X11,X12,group,'rb','v^',[],'off');
set(h1,'LineWidth',2)
legend('Control','Alcoholic','Location','NW')

% Classify a grid of measurements on the same scale:

[X,Y] = meshgrid(linspace(min(X11),max(X11))*1.2,linspace(min(X12),max(X12))*1.2);
X = X(:); Y = Y(:);
[C,err,P,logp,coeff] = classify([X Y],[X11 X12],group,'quadratic');

% Visualize the classification:

hold on;
gscatter(X,Y,C,'rb','.',1,'off');
K = coeff(1,2).const;
L = coeff(1,2).linear; 
Q = coeff(1,2).quadratic;
% Function to compute K + L*v + v'*Q*v for multiple vectors
% v=[x;y]. Accepts x and y as scalars or column vectors.
f = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);

h2 = ezplot(f,[min(X11) max(X11) min(X12) max(X12)]*1.2);
set(h2,'Color','m','LineWidth',2)
title('{\bf Classification with EEG Data}')

end

%profile viewer;