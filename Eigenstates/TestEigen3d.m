clc;
clear;
close all;
a=10; hbar=1; m=1;dx=2/13+10^-4;
x=[-a:dx:a]'; y=x;z=x;[XX,YY,ZZ]=meshgrid(z,y,x);
V=-5/sqrt(XX.^2+YY.^2+ZZ.^2);
[W,D,flag]=Eigen3D(dx,m,hbar,V);
[W,R]=mgsog(W);
nx=length(x); ny=length(y);nz=length(z); 
%%
t=[0:1:48];
fig=figure(1); s=1;
pause(2)
for i=1:30;
    for k=1:length(t);
clf(fig)
colors=sign(reshape(W(:,i),nx-2,ny-2,nz-2));
isosurface(XX(2:nx-1,2:ny-1,2:nz-1),YY(2:nx-1,2:ny-1,2:nz-1),ZZ(2:nx-1,2:ny-1,2:nz-1),reshape(abs(W(:,i)).^2,nx-2,ny-2,nz-2),colors)
view(-100+5.42*t(k),15)
camlight('right')
colormap([0 0 1; 1 0.15 0])
axis equal
axis off
set(gcf,'InvertHardCopy','off','Color','white')
F(s)=getframe(gcf);
s=s+1;
    end
end
%%
video=VideoWriter('Eigen3DHy2');
video.Quality=100;
video.FrameRate=24;
open(video);
writeVideo(video,F);
close(video);
%%
tn=[0:0.01:5]; s=1;
fig2=figure(2);
pause(2)
for k=1:length(tn);
Wt=1/sqrt(2)*W(:,10)*exp(-j*D(10,10)*tn(k)/hbar)+j/sqrt(2)*W(:,27)*exp(-j*D(27,27)*tn(k)/hbar);
colors=sign(reshape(real(Wt),nx-2,ny-2,nz-2));
isosurface(XX(2:nx-1,2:ny-1,2:nz-1),YY(2:nx-1,2:ny-1,2:nz-1),ZZ(2:nx-1,2:ny-1,2:nz-1),reshape(abs(Wt).^2,nx-2,ny-2,nz-2),colors)
view(-100+45*tn(k),70)
camlight('right')
colormap([0 0 1; 1 0.15 0])
axis([-a/1.6 a/1.6 -a/1.6 a/1.6 -a/2 a/2])
axis off
set(gcf,'InvertHardCopy','off','Color','white')
Fa(s)=getframe(gcf);
clf(fig2)
s=s+1;
end
%%
video=VideoWriter('Eigen3DHAnim3');
video.Quality=100;
video.FrameRate=24;
open(video);
writeVideo(video,Fa);
close(video);