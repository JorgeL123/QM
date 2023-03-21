clc;
clear;
close all;
%%
dx=0.04; dt=0.01; hbar=1;m=10; a=10;
x=[-1.05*a:dx:1.05*a]; y=[-1.05*a:dx:1.05*a]; [XX,YY]=meshgrid(x,y); t=[0:dt:140];c=10; d=1;
P=1/d*(2/pi)^(1/2)*exp((-(XX).^2-(YY).^2)/d^2).*exp(j*c*YY); n=1;s=1; ang=[0:pi/40:2*pi]; z1=0.85; z2=0.6; %z=0.7;
for i=1:10:length(t);
    T(s)=t(i);
    s=s+1;
end
%Define potential
r=hbar*c/(2*m)*dt/dx; la=a/5; lb=a/7; lc=a/100;
xr=[-lc lc lc -lc]; yr1=[la la a a]; yr2=[-lb -lb lb lb]; yr3=[-a -a -la -la];
%%
 ang=[0:pi/80:2*pi];
Y=z1*cos(ang)'; Y=[(Y(1)+Y(end))/2;Y(2:end-1);(Y(1)+Y(end))/2];
Z=z2*sin(ang)'; Z=[(Z(1)+Z(end))/2;Z(2:end-1);(Z(1)+Z(end))/2];
Y=interp(Y,2); Z=interp(Z,2);
plot(Y*a,Z*a,'-b')
axis([-1.05*a 1.05*a -0.4*a 0.55*a 0 1])
f=(Y+j*Z)*a;
N=length(f);
%%
J=zeros(length(x));
in=inpolygon(XX,YY,Y*a,Z*a);
xf=XX(in); yf=YY(in);
for i=1:length(xf);
        J=(XX==xf(i)).*(YY==yf(i))+J;
end
imagesc(x,y,J)
hold on
plot(Y*a,Z*a,'-k','LineWidth',2)
set(gca,'YDir','normal')
hold off
J=abs(J-1);
%%
V=90*J;
[W,D,flag]=Eigen2D(dx,m,hbar,V);
%%
J=2*J-1;
J(J==1)=NaN; J(J==-1)=1;
nx=length(x); ny=length(y);
imagesc(x,y,J)
hold on
plot(Y*a,Z*a,'-k','LineWidth',2)
axis([-a a -0.5*a 0.5*a 0 1])
set(gca,'YDir','normal')
%%
W=-W;
%%
fh = figure('Menu','none','ToolBar','none'); 
ah = axes('Units','Normalize','Position',[0 0 1 1])
pause(1.5)
%%
base = [0 1 0; 0 1 0;0 1 0;0 1 0;0 1 0;0 1 0;0 1 0; 0 0 1;0 0 1;0 0 1; 1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0;1 0 0];
customMap = abs(interp1(linspace(length(x),0,size(base,1)), base, fliplr(0:length(x)), 'pchip'));
for i=1:3
    for s=1:length(x)
        if customMap(s,i)>1;
            customMap(s,i)=1;
        end
    end
end
video=VideoWriter('PondEi');
video.Quality=100;
video.FrameRate=40;
open(video);
for i=1:length(W(1,:))-1;
 for j=1:24;
pcolor(XX(2:nx-1,2:ny-1),YY(2:nx-1,2:ny-1),real(reshape(W(:,i),nx-2,ny-2)).*J(2:nx-1,2:ny-1))
set(gca,'YDir','normal')
colormap(customMap)
shading interp
axis([-1.2*a 1.2*a -0.75*a 0.75*a 0 1])
xlabel('x')
ylabel('y')
zlabel('Ψ')
hold on
%legend('V(x)','|Ψ(x,t)|','Re\{Ψ(x,t)\}')`
plot(Y*a,Z*a,'-k','LineWidth',6.2)
set(gcf,'InvertHardCopy','off','Color','black');
clim([-max(max(W)) max(max(W))])
text(0.95*a,-0.45*a,strcat(num2str(i)),'interpreter','latex','FontSize',40,'color','white')
axis off
print(gcf,'foo.png','-dpng','-r300');
A=imread('foo.png');
writeVideo(video,A);
 end
 hold off
 for k=0:24;
     pcolor(XX(2:nx-1,2:ny-1),YY(2:nx-1,2:ny-1),((24-k)*real(reshape(W(:,i),nx-2,ny-2))+k*real(reshape(W(:,i+1),nx-2,ny-2)))/(24).*J(2:nx-1,2:ny-1))
set(gca,'YDir','normal')
colormap(customMap)
shading interp
axis([-1.2*a 1.2*a -0.75*a 0.75*a 0 1])
xlabel('x')
ylabel('y')
zlabel('Ψ')
hold on
plot(Y*a,Z*a,'-k','LineWidth',6.2)
set(gcf,'InvertHardCopy','off','Color','black');
clim([-max(max(W)) max(max(W))])
axis off
print(gcf,'foo.png','-dpng','-r300');
A=imread('foo.png');
writeVideo(video,A);
 end
end
close(video);
%%
pcolor(XX(2:nx-1,2:ny-1),YY(2:nx-1,2:ny-1),real(reshape(W(:,80),nx-2,ny-2)).*J(2:nx-1,2:ny-1))
set(gca,'YDir','normal')
colormap(customMap)
shading interp
axis([-1.2*a 1.2*a -0.75*a 0.75*a 0 1])
xlabel('x')
ylabel('y')
zlabel('Ψ')
hold on
%legend('V(x)','|Ψ(x,t)|','Re\{Ψ(x,t)\}')`
plot(Y*a,Z*a,'-k','LineWidth',6.2)
set(gcf,'InvertHardCopy','off','Color','black');
text(0.95*a,-0.45*a,strcat(num2str(80)),'interpreter','latex','FontSize',40,'color','white')
clim([-max(max(W)) max(max(W))])
axis off
print(gcf,'foo2.png','-dpng','-r300');