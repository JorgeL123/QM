function [W,D,A,flag]=Eigen1D(n,dx,m, hbar,V);
V=V(2:n-1,1);
n=n-2;
d1=ones(n,1)*(hbar^2)/(m*dx^2)+V;
d2=-ones(n,1)*(hbar^2)/(2*m*dx^2);
A=spdiags([d2 d1 d2],-1:1,n,n);
[W,D,flag]=eigs(A,15,'smallestreal');
end
