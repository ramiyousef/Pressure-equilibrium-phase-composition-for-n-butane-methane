function fp=dertermproject(fun,x,h,zi,Ki,N)
% Calculates the derivate of a function numerically based
% on the 5-pt Stencil
% (c) 2011 Phillip Servio
fp=(-feval(fun,x+2*h,zi,Ki,N)+8*feval(fun,x+h,zi,Ki,N)-8*feval(fun,x-h,zi,Ki,N)+feval(fun,x-2*h,zi,Ki,N))/(12*h);
end