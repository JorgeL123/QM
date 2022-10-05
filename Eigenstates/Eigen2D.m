function [W,D,flag]=Eigen2D(dx,m,hbar,V);
[nx,ny]=size(V);
V=V(2:nx-1,2:ny-1);
nx=nx-2; ny=ny-2;
V=reshape(V,nx*ny,1);
c=(hbar^2)/(2*m*dx^2);
da=-ones(nx*ny,2);
da(1:nx:nx*ny,:)=0;
diagonals=[4*ones(nx*ny,1)+V/c,da,-ones(nx*ny,2)];
A=spdiags(diagonals,[0 -1 1 -nx nx],nx*ny,nx*ny);
[W,D,flag]=eigs(A,15,'smallestreal');
D=c*D;
end