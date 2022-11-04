clc;
clear;
close all;
dx=0.02; dt=0.01; hbar=1;m=125; a=10;kappa=0.1;
x=[-2*a:dx:2*a]; t=[0:dt:400];c=6;
P=(2/pi)^(1/4)*exp(-(x+a/2).^2).*exp(j*c*x); T=[1:0.75:500]; n=1; V0=0.3;
%Define potential
for i=1:length(x);
if abs(x(i))<=a/5&&abs(x(i))>=a/15;
V(i)=V0;
elseif abs(x(i))>=a/5;
V(i)=0;
elseif abs(x(i))<a/15;
V(i)=0;
end
end
V=0*x;
%V=(x.^4/35-x.^2+4)/17-0.3;
d=0.009;
V=0.1*(1/d*exp(-(x).^2/d^2))+1/2*1/25*x.^2;
r=hbar*c/(2*m)*dt/dx;
%%
figure(1);
pause(1.5)
for k=1:length(t);
 % Draw Wave Function for certain values of time
if t(k)==T(n);
plot(x,V,'-k','LineWidth',1.1)
hold on
plot(x,abs(P),'-b','LineWidth',1.5)
hold on
plot(x,real(P),'-r','LineWidth',1.5)
axis([-a a -1.1 1.1])
xlabel('$x$','interpreter','latex')
ylabel('$\Psi$','interpreter','latex')
legend('$V(x)$','$|\Psi(x,t)|$','$Re\{\Psi(x,t)\}$','interpreter','latex')
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])
F(n)=getframe(gcf);
drawnow
n=n+1;
hold off
end
% RK4
K1=NSchroE(P,V,dx,hbar,m,kappa);
K2=NSchroE(P+dt*K1/2,V,dx,hbar,m,kappa);
K3=NSchroE(P+dt*K2/2,V,dx,hbar,m,kappa);
K4=NSchroE(P+dt*K3,V,dx,hbar,m,kappa);
A=P(end); B=P(end-1);
C=P(1); D=P(2);
P=P+dt/6*(K1+2*K2+2*K3+K4);
% Mur Boundary condition: not very clean for big wave packets.
P(end)=B+(r-1)/(r+1)*(P(end-1)-A);
P(1)=D+(r-1)/(r+1)*(P(2)-C);
end
%%
video=VideoWriter('NLSE8');
video.Quality=100;
video.FrameRate=24;
open(video);
writeVideo(video,F);
close(video);