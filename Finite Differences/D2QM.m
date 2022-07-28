clc;
clear;
close all;
dx=0.1; dt=0.01; hbar=1;m=10; a=10;
x=[-1.5*a:dx:1.5*a]; y=[-1.5*a:dx:1.5*a]; [XX,YY]=meshgrid(x,y); t=[0:dt:50];c=2;
P=(2/pi)^(1/2)*exp(-(XX).^2-(YY+a/3).^2).*exp(j*c*YY/2+j*c*XX/2); T=[0:0.25:50]; n=1; alfa=[0:pi/40:2*pi];
%Define potential
V=(XX.^2+YY.^2)/20;
r=hbar*c/(2*m)*dt/dx;
%%
figure(1);
pause(1.5)
for k=1:length(t);
 % Draw Wave Function for certain values of time
if t(k)==T(n);
%plot3(reshape(XX,length(x).^2,1),reshape(YY,length(x).^2,1),reshape(V,length(x).^2,1),'-k','LineWidth',.001)
plot3(-y,zeros(length(x),1),(y.^2)/20,'-k','LineWidth',1.3)
hold on
hold on
plot3(zeros(length(x),1),x,(x.^2)/20,'-k','LineWidth',1.3)
hold on
plot3(sqrt(20)*cos(alfa),sqrt(20)*sin(alfa),ones(length(alfa),1),'-k','LineWidth',1.3)
surf(XX,YY,abs(P),'FaceAlpha',0.9)
colormap jet
shading interp
axis([-a a -a a 0 1])
xlabel('x')
ylabel('y')
zlabel('Ψ')
%legend('V(x)','|Ψ(x,t)|','Re\{Ψ(x,t)\}')
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])
set(gca,'Ztick',[])
F(n)=getframe(gcf);
n=n+1;
hold off
end
% RK4
%P(1,:)=0; P(end,:)=0;P(:,1)=0; P(:,end)=0;
K1=SchroE2(P,V,dx,hbar,m);
K2=SchroE2(P+dt*K1/2,V,dx,hbar,m);
K3=SchroE2(P+dt*K2/2,V,dx,hbar,m);
K4=SchroE2(P+dt*K3,V,dx,hbar,m);
%A=P(end); B=P(end-1);
%C=P(1); D=P(2);
P=P+dt/6*(K1+2*K2+2*K3+K4);
% Mur Boundary condition: not very clean for big wave packets.
%P(end)=B+(r-1)/(r+1)*(P(end-1)-A);
%P(1)=D+(r-1)/(r+1)*(P(2)-C);
end
%%
video=VideoWriter('QM2D2');
video.Quality=100;
video.FrameRate=24;
open(video);
writeVideo(video,F);
close(video);
