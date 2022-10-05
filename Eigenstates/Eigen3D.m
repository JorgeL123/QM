function [W,D,flag]=Eigen3D(dx,m,hbar,V);
[nx,ny,nz]=size(V);
V=V(2:nx-1,2:ny-1,2:nz-1);
nx=nx-2; ny=ny-2; nz=nz-2;
V=reshape(V,nx*ny*nz,1);
c=(hbar^2)/(2*m*dx^2);
da=-ones(nx*ny*nz,2);
da(1:nx:nx*ny*nz,:)=0;
db=-ones(nx*ny*nz,2);
db(1:nx*ny:nx*ny*nz,:)=0;
diagonals=[6*ones(nx*ny*nz,1)+V/c,da,db,-ones(nx*ny*nz,2)];
A=spdiags(diagonals,[0 -1 1 -nx nx -nx*ny nx*ny],nx*ny*nz,nx*ny*nz);
[W,D,flag]=eigs(A,30,'smallestreal');
D=c*D;
end