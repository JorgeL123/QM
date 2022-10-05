clc;
clear;
close all;
a=10; hbar=1; m=10;dx=0.071203445;
x=[-1.25*a:dx:1.25*a]'; y=x;[XX,YY]=meshgrid(x,y);
%V=0.001*(XX.^2+YY.^2);
V=log(sqrt(XX.^2+YY.^2)/a);
%V=((XX.^2+YY.^2).^2/70-XX.^2-YY.^2+26)/450-0.03;
[W,D,flag]=Eigen2D(dx,m,hbar,V);
%%
nx=length(x); ny=length(y);dt=1.5/10; tf=150/10;
t=[0:dt:tf];
figure(1);
pause(2)
n=0; s=1;%U=["Re\{\psi_{00}\}","Re\{\psi_{01}\}", "Re\{\psi_{10}\}", "Re\{\psi_{02}\}", "Re\{\psi_{20}\}", "Re\{\psi_{11}\}", "Re\{\psi_{03}\}", "Re\{\psi_{30}\}", "Re\{\psi_{21}\}", "Re\{\psi_{12}\}",...
   % "Re\{\psi_{04}\}", "Re\{\psi_{40}\}", "Re\{\psi_{31}\}", "Re\{\psi_{13}\}", "Re\{\psi_{22}\}"];
for g=1:15;
for k=1:length(t);
plot3(x,zeros(1,length(x)),log(x/a),'-k','LineWidth',1.3)
hold on
plot3(zeros(1,length(x)),x,log(x/a),'-k','LineWidth',1.3)
hold on
surf(XX(2:nx-1,2:ny-1),YY(2:nx-1,2:ny-1),real(reshape(W(:,g),nx-2,ny-2)*exp(-i*D(g,g)*t(k)/hbar)))
shading interp 
colormap jet
view(-45-(s)/5, 50)
%title([U(g)])
axis([-a/3 a/3 -a/3 a/3 -max(max((W(:,:)))) max(max(W(:,:)))])
set(gcf,'InvertHardCopy','off','Color','white');
axis off
F(s)=getframe(gcf);
hold off
s=s+1;
end
end
%%
video=VideoWriter('EigenHA2d');
video.Quality=100;
video.FrameRate=24;
open(video);
writeVideo(video,F);
close(video);