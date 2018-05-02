%  KDR line learch routine

function [New tr] = FoldedKDRline(X,n,Gy,sgx0,Q,L,R,dLR,ETA,EPS,index)

    function t=kdr1dim1(s)
        tmpR=R-s*dLR;
        [tmpR,D,V]=svd(tmpR,0);
        tmpB=kron(tmpR,L); 
        Z=X*tmpB;
        zzaa=sum(Z.*Z,2);
        zzab=Z*Z';
        D=repmat(zzaa,1,n);
        zz=abs(D + D' - 2*zzab);
        Kz=exp(-zz./sgx0);  
        mK=mean(Kz,2);
        rK=repmat(mK,1,n);
        Gz=Kz-rK-rK'+mean(mK,1)*ones(n,n);   
        Gz=Q*Gz*Q;
        Gz=(Gz+Gz')./2;
        %Gzi=inv(Gz+n*EPS.*eye(n));
        %t=sum(sum(Gy.*Gzi,1),2);
        t=trace((Gz+n*EPS.*eye(n))\Gy);
    end
 
    function t=kdr1dim2(s)
        tmpL=L-s*dLR;
        [tmpL,D,V]=svd(tmpL,0);
        tmpB=kron(R,tmpL); 
        Z=X*tmpB;
        zzaa=sum(Z.*Z,2);
        zzab=Z*Z';
        D=repmat(zzaa,1,n);
        zz=abs(D + D' - 2*zzab);
        Kz=exp(-zz./sgx0);  
        mK=mean(Kz,2);
        rK=repmat(mK,1,n);
        Gz=Kz-rK-rK'+mean(mK,1)*ones(n,n);     
        Gz=Q*Gz*Q;
        Gz=(Gz+Gz')./2;
        %Gzi=inv(Gz+n*EPS.*eye(n));
        %t=sum(sum(Gy.*Gzi,1),2);
        t=trace((Gz+n*EPS.*eye(n))\Gy);
   end

if index==1
        
[s tr] = fminbnd(@kdr1dim1,0,ETA,optimset('Display', 'off', 'MaxIter',100));
New = R-s*dLR;
[New,D,V]=svd(New,0);

else if index==2

[s tr] = fminbnd(@kdr1dim2,0,ETA,optimset('Display', 'off', 'MaxIter',100));
New = L-s*dLR;  
[New,D,V]=svd(New,0);

    end
end



end