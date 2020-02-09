function [w] = GCL_FG(S,a,i,j,g,h,n,m,index,RB)
xi=S.x(i,:); % central point
xj=S.x(j,:);
ui=S.u(i,:); % directon vector
uj=S.u(j,:);
hj=S.h(j); % height
rj=S.r(j);
a1 = hj*uj(3)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2))+sqrt(ui(2)^2+ui(3)^2)*(-xi(1)+xj(1))+ui(1)*(xi(3)-xj(3));
a2 = uj(2)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2));
a3 = ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2);
a4 =  hj*(-uj(2)*ui(3)+ui(2)*uj(3)*(ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)))/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))+ui(1)*ui(2)*(-xi(1)+xj(1))/sqrt(ui(2)^2+ui(3)^2)+ui(3)*(-xi(2)+xj(2))/sqrt(ui(2)^2+ui(3)^2)+ui(2)*(-xi(3)+xj(3));
a5 = ui(2)*(-uj(1)+ui(1)*sqrt(uj(2)^2+uj(3)^2)/sqrt(ui(2)^2+ui(3)^2));
a6 = (ui(1)*uj(1)*ui(2)*uj(2)+ui(3)*uj(3)+ui(2)*uj(2)*sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2));
b1 = hj*uj(3)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2))+sqrt(ui(2)^2+ui(3)^2)*(-xi(1)+xj(1))+ui(1)*(xi(3)-xj(3));
b2 = uj(2)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2));
b3 = ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2);
c1 = hj*(ui(2)*uj(2)+ui(3)*uj(3)*(ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)))/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))+ui(1)*ui(3)*(-xi(1)+xj(1))/sqrt(ui(2)^2+ui(3)^2)+ui(2)*(xi(2)-xj(2))/sqrt(ui(2)^2+ui(3)^2)+ui(3)*(-xi(3)+xj(3));
c2 = ui(3)*(-uj(1)+ui(1)*sqrt(uj(2)^2+uj(3)^2)/sqrt(ui(2)^2+ui(3)^2));
c3 = (ui(1)*uj(1)*uj(2)*ui(3)-ui(2)*uj(3)+uj(2)*ui(3)*sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2));
%Qi = @(x,y) sqrt( ( hj*uj(3)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2))+sqrt(ui(2)^2+ui(3)^2)*(-xi(1)+xj(1))+ui(1)*(xi(3)-xj(3))+uj(2)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2)).*x.*sin(y)+(ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)).*x.*cos(y)).^2 + ( hj*(-uj(2)*ui(3)+ui(2)*uj(3)*(ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)))/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))+ui(1)*ui(2)*(-xi(1)+xj(1))/sqrt(ui(2)^2+ui(3)^2)+ui(3)*(-xi(2)+xj(2))/sqrt(ui(2)^2+ui(3)^2)+ui(2)*(-xi(3)+xj(3))+ui(2)*(-uj(1)+ui(1)*sqrt(uj(2)^2+uj(3)^2)/sqrt(ui(2)^2+ui(3)^2)).*x.*cos(y)+(ui(1)*uj(1)*ui(2)*uj(2)+ui(3)*uj(3)+ui(2)*uj(2)*sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)).*x.*sin(y)/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))).^2 );
%Wi = @(x,y) acos( (hj*uj(3)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2))+sqrt(ui(2)^2+ui(3)^2)*(-xi(1)+xj(1))+ui(1)*(xi(3)-xj(3))+uj(2)*(-ui(1)+uj(1)*sqrt(ui(2)^2+ui(3)^2)/sqrt(uj(2)^2+uj(3)^2)).*x.*sin(y)+(ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)).*x.*cos(y))./Qi(x,y) );
%Ei = @(x,y) hj*(ui(2)*uj(2)+ui(3)*uj(3)*(ui(1)*uj(1)+sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)))/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2))+ui(1)*ui(3)*(-xi(1)+xj(1))/sqrt(ui(2)^2+ui(3)^2)+ui(2)*(xi(2)-xj(2))/sqrt(ui(2)^2+ui(3)^2)+ui(3)*(-xi(3)+xj(3))+ui(3)*(-uj(1)+ui(1)*sqrt(uj(2)^2+uj(3)^2)/sqrt(ui(2)^2+ui(3)^2)).*x.*cos(y)+(ui(1)*uj(1)*uj(2)*ui(3)-ui(2)*uj(3)+uj(2)*ui(3)*sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2)).*x.*sin(y)/(sqrt(ui(2)^2+ui(3)^2)*sqrt(uj(2)^2+uj(3)^2));
Qi = @(x,y) sqrt((a1+a2.*x.*sin(y)+a3.*x.*cos(y)).^2+(a4+a5.*x.*cos(y)+a6.*x.*sin(y)).^2);
Wi = @(x,y) acos( (b1+b2.*x.*sin(y)+b3.*x.*cos(y))./Qi(x,y) );
Ei = @(x,y) c1+c2.*x.*cos(y)+c3.*x.*sin(y);
if(index==0) %F
    if(n==0)   
        %f = @(x,y)  x.*besselj(0,GCL_rootBessel(1,m).*x).*besselj(g,GCL_rootBessel(g+1,h).*Qi(x,y)).*sin(g*Wi(x,y)).*exp(-GCL_rootBessel(g+1,h).*Ei(x,y));
        f = @(x,y) x.*besselj(0,RB(1,m).*x).*besselj(g,RB(g+1,h).*Qi(x,y)).*sin(g*Wi(x,y)).*exp(-RB(g+1,h).*Ei(x,y));
        w = integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(1,RB(1,m))^2)*exp(-RB(1,m)*hj));
    else
        f = @(x,y) x.*besselj(n,RB(n+1,m).*x).*cos(n.*y).*besselj(g,RB(g+1,h).*Qi(x,y)).*sin(g.*Wi(x,y)).*exp(-RB(g+1,h).*Ei(x,y));
        w = 2*integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(n+1,RB(n+1,m))^2)*exp(-RB(1,m)*hj));
    end
end
if(index==1) %G
    if(n==0)
        f = @(x,y) x.*besselj(0,RB(1,m).*x).*besselj(g,RB(g+1,h).*Qi(x,y)).*cos(g.*Wi(x,y)).*exp(-RB(g+1,h).*Ei(x,y));
        w = integral2(f,0,rj,0,2*pi)/(pi*(a^2)*besselj(1,RB(1,m))^2*exp(-RB(1,m)*hj));
    else
        f = @(x,y) x.*besselj(n,RB(n+1,m).*x).*cos(n.*y).*besselj(g,RB(g+1,h).*Qi(x,y)).*cos(g.*Wi(x,y)).*exp(-RB(g+1,h).*Ei(x,y));
        w = 2*integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(n+1,RB(n+1,m))^2)*exp(-RB(1,m)*hj));
    end
end
if(index==2) %F'
    if(n==0)
        w = 0;
    else
        f = @(x,y) x.*besselj(n,RB(n+1,m).*x).*sin(n.*y).*besselj(g,RB(g+1,h).*Qi(x,y)).*sin(g.*Wi(x,y)).*exp(-RB(g+1,h).*Ei(x,y));
        w = 2*integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(n+1,RB(n+1,m))^2)*exp(-RB(1,m)*hj));
    end
end
if(index==3) %G'
    if(n==0)
        w = 0;
    else
        f = @(x,y) x.*besselj(n,RB(n+1,m).*x).*sin(n.*y).*besselj(g,RB(g+1,h).*Qi(x,y)).*cos(g.*Wi(x,y)).*exp(-RB(g+1,h).*Ei(x,y));
        w = 2*integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(n+1,RB(n+1,m))^2)*exp(-RB(1,m)*hj));
    end
end