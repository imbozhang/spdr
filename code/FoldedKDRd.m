function dNew= FoldedKDRd (X,B,L,R,n,sgx,EPS,Gy,pr,dr,pl,dl,index)
  
    Z=X*B;                          
    zzaa=sum(Z.*Z,2);
    zzab=Z*Z';
    D=repmat(zzaa,1,n);
    zz=abs(D + D' - 2*zzab);
    Kz=exp(-zz./sgx);
    mK=mean(Kz,2);
    rK=repmat(mK,1,n);
    Gz=Kz-rK-rK'+mean(mK,1)*ones(n,n); 
    Gze=Gz+EPS.*n*eye(n);
    GziGyGzi=Gze\Gy/Gze;  
    
if index==1
            
        Lsum=sum(L,2); % pl*1
            
        for a=1:pr
        Xaa=X(:,(a-1)*pl+1:a*pl);          
        Xa=repmat(Xaa*Lsum,1,n);
        XX=Xa-Xa';
        
            for b=1:dr
            Zsum=sum(Z(:,(b-1)*dl+1:b*dl),2);
            Zb=repmat(Zsum,1,n);
            ZZ=Zb-Zb';
            tt=XX.*ZZ.*Kz;
            mK=mean(tt,2);
            rK=repmat(mK,1,n);
            dKR=tt-rK-rK'+mean(mK,1)*ones(n,n);        
            dR(a,b)=trace(GziGyGzi*dKR);
            end
        end
        
         dNew=dR/norm(dR);
    
    else if index==2
    
        Rsum=sum(R,2); % pr*1
    
        for a=1:pl
        Xaa=X(:,a:pl:length(X(1,:)));
        Xa=repmat(Xaa*Rsum,1,n);
        XX=Xa-Xa';
            for b=1:dl
            if dr==1
            Zb=repmat(Z(:,b),1,n); 
                else if dr==2
            Zb=repmat(Z(:,b)+Z(:,dl+b),1,n); 
                       else if dr==3
            Zb=repmat(Z(:,b)+Z(:,dl+b)+Z(:,2*dl+b),1,n); 
                           end
                        end
            end
            ZZ=Zb-Zb';
            tt=XX.*ZZ.*Kz;
            mK=mean(tt,2);
            rK=repmat(mK,1,n);
            dKL=tt-rK-rK'+mean(mK,1)*ones(n,n);        
            dL(a,b)=trace(GziGyGzi*dKL);
            end
        end
        
    dNew=dL/norm(dL);
        end
end

            