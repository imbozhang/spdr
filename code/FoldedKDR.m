%  Folded Kernel Dimension Reduction
%
%  Author: Bo Zhang 
%  Affiliation: North Carolina State University
%  Version: 0.1
%  Date: May 1, 2011
%
%  (C) Copyright:  Bo Zhang
%
%-----------------------------------------------
% FoldedKDR()
%
% Arguments
%  X:  explanatory variables (input data)
%  Y:  response variables (teaching data)
%
% Return value(s)
%  L:
%  R:
%  B:  
%  tr:
%-----------------------------------------------


function [L R B tr]=FoldedKDR(X,Y,pl,dl,pr,dr,sgx0,sgy0)

% EEG without pcr
% T=5000;                % Number of iterations in foldedKDR method (without pcr)
% ANL=200;               % Tuning parameter (without pcr)
% EPS=1e-004;            % Regularization parameter for matrix inversion.

% EEG with pcr
T=50;                    % Number of iterations in foldedKDR method (with pcr)
ANL=10;                  % Tuning parameter (with pcr)
EPS=1e-004;              % Regularization parameter for matrix inversion.

% example 1
% T=100;                 % Number of iterations in foldedKDR method 
% ANL=20;                % Tuning parameter (with pcr)
% EPS=1e-004;            % Regularization parameter for matrix inversion.

% binary 1
%T=50;                   % Number of iterations in foldedKDR method 
%ANL=10;                 % Tuning parameter (with pcr)
%EPS=1e-004;             % Regularization parameter for matrix inversion.


DISPtr=0;                % trace for optimization process
DISPplot=1;              % Graphical display for optimization process
ETA=10;                  % Range of the golden ratio search
TOL=1e-008;              % Tolerance for stopping the optimization

L0=randn(pl,dl);         % L is initialized
R0=randn(pr,dr);         % R is initialized

[L,D,V]=svd(L0,0);       % Normalization
[R,D,V]=svd(R0,0);       % Normalization
B=kron(R,L);             % B is the direct product of beta and alpha

n=length(Y(:,1));
unit=ones(n,n);
I=eye(n);
Q=I-unit./n;

yyaa=sum(Y.*Y,2);
yyab=Y*Y';
D=repmat(yyaa,1,n);
yy=abs(D + D' - 2*yyab);

Ky=exp(-yy./sgy0);
mK=mean(Ky,2);
rK=repmat(mK,1,n);
Gy=Ky-rK-rK'+mean(mK,1)*unit;  
Gy=Q*Gy*Q;
Gy=(Gy+Gy')./2;


Z=X*B;
zzaa=sum(Z.*Z,2);
zzab=Z*Z';
D=repmat(zzaa,1,n);
zz=abs(D + D' - 2*zzab);
Kz=exp(-zz./sgx0);
mK=mean(Kz,2);
rK=repmat(mK,1,n);
Gz=Kz-rK-rK'+mean(mK,1)*unit; 
Gz=Q*Gz*Q;
Gz=(Gz+Gz')./2;

tr0=trace((Gz+n*EPS.*I)\Gy);

RGN=5;

if DISPplot
   fxy = figure;
   axis([-RGN RGN -RGN RGN]);
   axis square;
   KDRplot(fxy,X,Y,B);
   drawnow;
end

mintr=1e+008;

%% Main loop for minimization
for h=1:T
    sgx=sgx0+(ANL-1)*sgx0*(T-h)/T;
    sgy=sgy0+(ANL-1)*sgy0*(T-h)/T;
%% recompute Gy    
    Ky=exp(-yy./sgy);
    mK=mean(Ky,2);
    rK=repmat(mK,1,n);
    Gy=Ky-rK-rK'+mean(mK,1)*unit;  
    Gy=Q*Gy*Q;
    Gy=(Gy+Gy')./2;
       
%% Update L by Steepest Descent with Line Search   
    dL=FoldedKDRd(X,B,L,R,n,sgx,EPS,Gy,pr,dr,pl,dl,2);
    if norm(dL)< TOL 
        break
    end
    [L tmptrL]=FoldedKDRline(X,n,Gy,sgx,Q,L,R,dL,ETA,EPS,2);
    B=kron(R,L);                   
 
%% Update R by Steepest Descent with Line Search   
    dR=FoldedKDRd(X,B,L,R,n,sgx,EPS,Gy,pr,dr,pl,dl,1);
    if norm(dR)< TOL 
        break
    end
    
    [R tmptrR]=FoldedKDRline(X,n,Gy,sgx0,Q,L,R,dR,ETA,EPS,1);
    B=kron(R,L);  
    
%% recompute the trace using original sgx and sgy
    
    aa=sum(Y.*Y,2);
    ab=Y*Y';
    D=repmat(aa,1,n);
    yy=max(D + D' - 2*ab, zeros(n,n));
    Gy=exp(-yy./sgy0);  
    Kyo=Q*Gy*Q;
    Kyo=(Kyo+Kyo')./2;  
    
    Z=X*B;
    nZ=Z./sqrt(sgx0);
    aa=sum(nZ.*nZ,2);
    ab=nZ*nZ';
    D=repmat(aa,1,n);
    ZZ=max(D + D' - 2*ab, zeros(n,n));
    Gz=exp(-ZZ);    
    Kz=Q*Gz*Q;
    Kz=(Kz+Kz')./2;     % Kz(i,j)=Q*exp(-||B^TX(i)-B^TX(j)||^2)*Q

    mz=inv(Kz+EPS*n.*I);
    tr=sum(sum(Kyo.*mz,1),2);
    
    if DISPtr
    fprintf('[%d]  tr=%.10f\n', h, tr);   
    end
    
    if tr < mintr
    mintr=tr;
    minL =L;
    minR =R;
    minB =B;
    end
    
    if DISPplot     
       KDRplot(fxy,X,Y,B);
       drawnow;
    end

end

L=minL;
R=minR;
B=minB;
tr=mintr;
