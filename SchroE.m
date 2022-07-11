function dP=SchroE(P,V,dx,hbar,m);
%dP(1)=j*hbar/(2*m)*(P(2)-2*P(1))/dx^2-j*V(1)*P(1)/hbar;
%dP(length(P))=j*hbar/(2*m)*(-2*P(length(P))+P(length(P)-1))/dx^2-j*V(length(P))*P(length(P))/hbar;
dP(1)=0; dP(length(P))=0;
%
dP(2)=j*hbar/(2*m)*(P(1)-2*P(2)+P(3))/(dx^2)-j*V(2)*P(2)/hbar;
dP(length(P)-1)=j*hbar/(2*m)*(P(length(P))-2*P(length(P)-1)+P(length(P)-2))/(dx^2)-j*V(length(P)-1)*P(length(P)-1)/hbar;
%
for i=3:length(P)-2;
    dP(i)=j*hbar/(2*m)*(-P(i+2)+16*P(i+1)-30*P(i)+16*P(i-1)-P(i-2))/(12*dx^2)-j*V(i)*P(i)/hbar;
end

end