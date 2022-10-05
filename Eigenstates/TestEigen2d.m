clc;
clear;
close all;
a=10; hbar=1; m=3;dx=0.1024162801;
x=[-a:dx:a]'; y=x;[XX,YY]=meshgrid(y,x);
%V=5*(XX.^2+YY.^2);
%V=-8./sqrt(XX.^2+YY.^2);
%V=((XX.^2+YY.^2).^2/25-XX.^2-YY.^2+16)/400-0.03;
V=((XX.^2+YY.^2).^2/30-XX.^2-YY.^2+3)/150+0.03;
[W,D,flag]=Eigen2D(dx,m,hbar,V);
%%
nx=length(x); ny=length(y);dt=2; tf=200;
t=[0:dt:tf];
figure(1);
pause(2)
n=0; s=1;
for g=1:15;
for k=1:length(t);
plot3(x,zeros(1,length(x)),(((x.^2).^2/30-x.^2+3)/150+0.01)/1.25,'-k','LineWidth',1.3)
hold on
plot3(zeros(1,length(x)),x,(((x.^2).^2/30-x.^2+3)/150+0.01)/1.25,'-k','LineWidth',1.3)
hold on
surf(XX(2:nx-1,2:ny-1),YY(2:nx-1,2:ny-1),real(reshape(W(:,g),nx-2,ny-2)*exp(-i*D(g,g)*t(k)/hbar)))
shading interp 
colormap jet
view(-45-(s)/5, 50)
title(['Re\{\psi_{',num2str(g),'}\}'])
axis([-a/1.5 a/1.5 -a/1.5 a/1.5 -max(max((W(:,:))))-0.006 max(max((W(:,:))))+0.006])
set(gcf,'InvertHardCopy','off','Color','white');
axis off
F(s)=getframe(gcf);
hold off
s=s+1;
end
end
%%
video=VideoWriter('EigenDuffing2df');
video.Quality=100;
video.FrameRate=24;
open(video);
writeVideo(video,F);
close(video);