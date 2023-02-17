clc;
clear;
close all;
load Data001
%%
dx=0.045; dt=0.01; hbar=1;m=10; a=10;
x=[-0.9*a:dx:0.9*a]; y=[-0.9*a:dx:0.9*a]; [XX,YY]=meshgrid(x,y); t=[0:dt:140];c=10; d=1;
P=1/d*(2/pi)^(1/2)*exp((-(XX).^2-(YY).^2)/d^2).*exp(j*c*YY); n=1;s=1; ang=[0:pi/40:2*pi]; z1=0.85; z2=0.6; %z=0.7;
for i=1:10:length(t);
    T(s)=t(i);
    s=s+1;
end
%Define potential
r=hbar*c/(2*m)*dt/dx; la=a/5; lb=a/7; lc=a/100;
xr=[-lc lc lc -lc]; yr1=[la la a a]; yr2=[-lb -lb lb lb]; yr3=[-a -a -la -la];
%%
Data=Data001; ang=[0:pi/10:2*pi];
Data(:,1)=smooth(Data(:,1));Data(:,2)=smooth(Data(:,2));
Y=[Data(4:end,1);-flip(Data(4:end,1))]; Y=[(Y(1)+Y(end))/2;Y(2:end-1);(Y(1)+Y(end))/2];
Z=[Data(4:end,2);flip(Data(4:end,2))]; Z=[(Z(1)+Z(end))/2;Z(2:end-1);(Z(1)+Z(end))/2];
Y=interp(Y,2); Z=interp(Z,2);
plot(Y*a,Z*a,'-b')
axis([-1.05*a 1.05*a -0.4*a 0.55*a 0 1])
f=(Y+j*Z)*a;
N=length(f);
%%
for k=1:N;
    s=0;
    for i=1:N;
        s=s+f(i)*exp(-2*pi*(k-1)*(i-1)*j/N)/N;
    end
    X(k)=s;
end
%%
fh = figure('Menu','none','ToolBar','none'); 
ah = axes('Units','Normalize','Position',[0 0 1 1])
pause(2);
video=VideoWriter('BatmanA');
video.Quality=100;
video.FrameRate=60;
open(video);
for k=1:N;
    s=0;
    for i=1:N;
        sa=s;
        s=s+X(i)*exp(2*pi*(k-1)*(i-1)*j/N);
        if round(i/2)==i/2;
        plot([real(sa);real(s)],[imag(sa);imag(s)],'-r','LineWidth',1.3)
        
        else
        plot([real(sa);real(s)],[imag(sa);imag(s)],'-y','LineWidth',1.3)
        end
        hold on
    end
    plot(0,0,'.g','MarkerSize',15)
    G(k)=real(s);
    F(k)=imag(s);
    plot(G(1:k),F(1:k),'-k','LineWidth',2)
    
       axis([-1.05*a 1.05*a -0.5*a 0.5*a 0 1])
        %set(gca,'Ytick',[]) 
%set(gca,'Xtick',[])
 set(gcf,'InvertHardCopy','off','Color','blue');
 axis off
print(gcf,'foo.png','-dpng','-r300');
A=imread('foo.png');
writeVideo(video,A);
         hold off
end
%%
close(video);
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
%%
V=0*XX;
fh = figure('Menu','none','ToolBar','none'); 
ah = axes('Units','Normalize','Position',[0 0 1 1])
pause(1.5)
base = [0 0 1;0 0 1;0.24 0.3 0.98; 0.3010 0.7450 0.9330;0.1 0.95 0.4330; 0 1 0; 1 1 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0; 1 0 0;0.9 0 0;0.9 0 0;0.8 0 0;0.7 0 0]; %green, yellow, red
customMap = abs(interp1(linspace(length(x),0,size(base,1)), base, fliplr(0:length(x)), 'pchip'));
for i=1:3
    for s=1:length(x)
        if customMap(s,i)>1;
            customMap(s,i)=1;
        end
    end
end
video=VideoWriter('BatmanC');
video.Quality=100;
video.FrameRate=24;
open(video);
for k=1:length(t);
 % Draw Wave Function for certain values of time
 %P(abs(y)>la,abs(x)<lc)=0;
%P(abs(y)<lb,abs(x)<lc)=0;
if T(n)==t(k);
%imagesc(x,y,abs(P))
pcolor(XX,YY,abs(P))
set(gca,'YDir','normal')
colormap(customMap)
shading interp
axis([-1.05*a 1.05*a -0.5*a 0.5*a 0 1])
xlabel('x')
ylabel('y')
zlabel('Ψ')
hold on
%legend('V(x)','|Ψ(x,t)|','Re\{Ψ(x,t)\}')`
plot(Y*a,Z*a,'-k','LineWidth',2)
set(gcf,'InvertHardCopy','off','Color','blue');
axis off
print(gcf,'foo.png','-dpng','-r300');
A=imread('foo.png');
writeVideo(video,A);
drawnow
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
P=P.*J;
% Mur Boundary condition: not very clean for big wave packets.
%P(end)=B+(r-1)/(r+1)*(P(end-1)-A);
%P(1)=D+(r-1)/(r+1)*(P(2)-C);
end
%%
close(video);