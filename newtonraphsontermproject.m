function [xn,nrfail] = newtonraphsontermproject(xl,xr,z,k,N)
% Newton Raphson method for calculating root in the
% interval xrg=[xl xr]
% (c) 2012 Phillip Servio
%% Initializing
 yl=Palpha(xl,z,k,N);
 yr=Palpha(xr,z,k,N);
 xg=(xl*yr-xr*yl)/(yr-yl);
i=0;
nrfail=0;
check=1;
h=1e-4;
tol=1e-8;
%% Loop
while tol<check
    i=i+1;
    fp=dertermproject(@Palpha,xg,h,z,k,N);
    y=Palpha(xg,z,k,N);
    xn=xg-y/fp;
    if xn<0 || xn>1 || i>10 || fp==0
        nrfail=1;
%         error('newton raphson failed');
        break 
    end 
    if abs(xn) <= 1
        check=abs(xg-xn);
    else
    end
    xg=xn;
    if nrfail==1
        break
    end
end