function dP=SchroE2(P,V,dx,hbar,m);
dP=zeros(length(P(1,:)));
dP(1,:)=0; dP(end,:)=0;dP(:,1)=0; dP(:,end)=0;
%
for i=2:length(P(1,:))-1;
dP(2,i)=j*hbar/(2*m)*(P(1,i)-4*P(2,i)+P(2,i+1)+P(2,i-1))/(dx^2)-j*V(2,i)*P(2,i)/hbar;
dP(end-1,i)=j*hbar/(2*m)*(P(end,i)-4*P(end-1,i)+P(end-2,i)+P(end-1,i+1)+P(end-1,i-1))/(dx^2)-j*V(end-1,i)*P(end-1,i)/hbar;
dP(i,2)=j*hbar/(2*m)*(P(i,1)-4*P(i,2)+P(i+1,2)+P(i-1,2))/(dx^2)-j*V(i,2)*P(i,2)/hbar;
dP(i,end-1)=j*hbar/(2*m)*(P(i,end)-4*P(i,end-1)+P(i,end-2)+P(i+1,end-1)+P(i-1,end-1))/(dx^2)-j*V(i,end-1)*P(i,end-1)/hbar;
end
%
for i=3:length(P(1,:))-2;
    for k=3:length(P(1,:))-2;
    dP(i,k)=j*hbar/(2*m)*(-P(i+2,k)+16*P(i+1,k)-60*P(i,k)+16*P(i-1,k)-P(i-2,k)-P(i,k+2)+16*P(i,k+1)+16*P(i,k-1)-P(i,k-2))/(12*dx^2)-j*V(i,k)*P(i,k)/hbar;
    end
end

end