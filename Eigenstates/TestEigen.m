clc;
clear;
close all;
a=5; hbar=1; m=10;dx=0.03;
x=[-a:dx:a]';
V=0*(x).^2; V(end)=100; V(1)=100; V0=0.1;
for i=1:length(x);
if abs(x(i))<a/3 && abs(x(i))>a/8;
    V(i)=V0;
end
end
[W,D,flag]=Eigen1D(length(x),dx,m,hbar,V);
%%
figure(1)
plot(x,[0;W(:,1);0],'-b','LineWidth',1.3)
hold on
plot(x,[0;W(:,2);0],'-r','LineWidth',1.3)
hold on
plot(x,[0;W(:,3);0],'-c','LineWidth',1.3)
hold on
plot(x,[0;W(:,4);0],'-g','LineWidth',1.3)
hold off
axis([-a a -0.13 0.13])
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])
xlabel('x')
ylabel('\psi_{n}')
legend('\psi_{1}','\psi_{2}','\psi_{3}','\psi_{4}','best')
pause(1.6)
saveas(gcf,'FIW1.png')
%%
figure(2)
t=[0:0.8:480]; A=[0;W(:,1);0];B=[0;W(:,2);0];C=[0;W(:,3);0];E=[0;W(:,4);0];
pause(2)
for i=1:length(t);
 Psi=sqrt(1/2)*A*exp(-j*D(1,1)*t(i)/hbar)+sqrt(1/4)*B*exp(-j*D(2,2)*t(i)/hbar)+sqrt(1/8)*C*exp(-j*D(3,3)*t(i)/hbar)+sqrt(1/8)*E*exp(-j*D(4,4)*t(i)/hbar);
 plot(x,V,'-k','LineWidth',1.2)
 hold on
 plot(x,abs(Psi),'-b','LineWidth',1.3)
 hold on
 plot(x,real(Psi),'-r','LineWidth',1.3)
 legend('V(x)','|\Psi(x,t)|','Re{\{\Psi(x,t)\}}')
 hold off
 axis([-a a -0.13 0.13])
set(gca,'Ytick',[]) 
set(gca,'Xtick',[])
xlabel('x')
ylabel('\Psi')
F(i)=getframe(gcf);
end
%%
video=VideoWriter('EigenFIW');
video.Quality=100;
video.FrameRate=40;
open(video);
writeVideo(video,F);
close(video);