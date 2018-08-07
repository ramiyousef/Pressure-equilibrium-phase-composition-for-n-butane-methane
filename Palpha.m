function y=Palpha(alpha,z,K,N)
y=0;
for i=1:N
y=y+z(i)*(K(i)-1)/(1+alpha*(K(i)-1));
end