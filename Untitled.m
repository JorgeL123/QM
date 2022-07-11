clc;
clear;
close all;
dx=0.05; dt=0.01; hbar=1;m=125; a=10;
x=[-2*a:dx:2*a]; t=[0:dt:300];c=10;
P=(2/pi)^(1/4)*exp(-(x+a/2).^2).*exp(j*c*x); T=[1:1:500]; n=1; V0=0.3;
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
%V=(x.^4/25-x.^2+4)/15-0.3;
%V=x.^2/40-0.5;
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
xlabel('x')
ylabel('Ψ')
legend('V(x)','|Ψ(x,t)|','Re\{Ψ(x,t)\}')
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])
F(n)=getframe(gcf);
n=n+1;
hold off
end
% RK4
K1=SchroE(P,V,dx,hbar,m);
K2=SchroE(P+dt*K1/2,V,dx,hbar,m);
K3=SchroE(P+dt*K2/2,V,dx,hbar,m);
K4=SchroE(P+dt*K3,V,dx,hbar,m);
A=P(end); B=P(end-1);
C=P(1); D=P(2);
P=P+dt/6*(K1+2*K2+2*K3+K4);
% Mur Boundary condition: not very clean for big wave packets.
P(end)=B+(r-1)/(r+1)*(P(end-1)-A);
P(1)=D+(r-1)/(r+1)*(P(2)-C);
end
%%
video=VideoWriter('Thing2');
video.Quality=100;
video.FrameRate=20;
open(video);
writeVideo(video,F);
close(video);