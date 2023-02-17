clc;
clear;
close all;
dx=0.045; dt=0.01; hbar=1;m=10; a=10;
x=[-0.9*a:dx:0.9*a]; y=[-0.9*a:dx:0.9*a]; [XX,YY]=meshgrid(x,y); t=[0:dt:140];c=10; d=2;
P=1/d*(2/pi)^(1/2)*exp((-(XX).^2-(YY).^2)/d^2).*exp(j*c*XX); n=1;s=1; ang=[0:pi/40:2*pi]; z1=0.85; z2=0.6; %z=0.7;
for i=1:10:length(t);
    T(s)=t(i);
    s=s+1;
end
%Define potential
V=0*(XX.^2+YY.^2)/20;
G=(XX/(z1*a)).^2+(YY/(z2*a)).^2<1;
r=hbar*c/(2*m)*dt/dx; la=a/5; lb=a/7; lc=a/100;
xr=[-lc lc lc -lc]; yr1=[la la a a]; yr2=[-lb -lb lb lb]; yr3=[-a -a -la -la];
%%
fh = figure('Menu','none','ToolBar','none'); 
ah = axes('Units','Normalize','Position',[0 0 1 1])
pause(1.5)
base = [0 0 1;0 0 1;0.3010 0.7450 0.9330; 0 1 0; 1 1 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0;0.9 0 0;0.9 0 0;0.8 0 0;0.7 0 0]; %green, yellow, red
customMap = abs(interp1(linspace(length(x),0,size(base,1)), base, fliplr(0:length(x)), 'pchip'));
for i=1:3
    for s=1:length(x)
        if customMap(s,i)>1;
            customMap(s,i)=1;
        end
    end
end
video=VideoWriter('Pond2');
video.Quality=100;
video.FrameRate=24;
open(video);
for k=1:length(t);
 % Draw Wave Function for certain values of time
 %P(abs(y)>la,abs(x)<lc)=0;
%P(abs(y)<lb,abs(x)<lc)=0;
if T(n)==t(k);
imagesc(x,y,abs(P))
colormap(customMap)
shading interp
axis([-1.3*a 1.3*a -0.8*a 0.8*a 0 1])
xlabel('x')
ylabel('y')
zlabel('Ψ')
hold on
%legend('V(x)','|Ψ(x,t)|','Re\{Ψ(x,t)\}')`
plot(z1*a*cos(ang),z2*a*sin(ang),'-k','LineWidth',2)
set(gcf,'InvertHardCopy','off','Color','blue');
axis off
print(gcf,'foo.png','-dpng','-r300');
A=imread('foo.png');
writeVideo(video,A);
%drawnow
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
P=P.*G;
% Mur Boundary condition: not very clean for big wave packets.
%P(end)=B+(r-1)/(r+1)*(P(end-1)-A);
%P(1)=D+(r-1)/(r+1)*(P(2)-C);
end
%%
close(video);