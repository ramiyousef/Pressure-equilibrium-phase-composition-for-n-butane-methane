function [xn,nrfail]=bisection_termproject(Palpha,z,K)

xl = 0;
xr = 1;
yl=feval(Palpha, xl, z, K,2); 
yr=feval(Palpha, xr, z, K,2);
dy=abs(yr-yl);
nrfail=0;
check=1;
tol=1e-8;

while check > tol
    
    xmid=(xl+xr)/2;
    ymid=feval(Palpha, xmid, z, K,2);
    
    if sign(yl*ymid)==1
        
        xl=xmid;
        yl=ymid;
        
    else
        
        xr=xmid;
        yr=ymid;
        
    end
    
    dynew=abs(yr-yl);
    
    if dynew>dy
        nrfail=1;
        break
    end
    
    if abs(xr)<=1
        
        check=abs(xr-xl);
        
    else
        
        check=abs(1-(xl/xr));
    end
    
    xn=xr;
    
end

   