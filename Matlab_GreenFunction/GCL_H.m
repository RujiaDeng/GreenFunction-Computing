function [v] = GCL_H(S,yq,a,j,n,m,index,RB)
xj=S.x(j,:); % central point
uj=S.u(j,:);% directon vector
hj=S.h(j); % height
rj=S.r(j);%radius
Fai = @(x,y) 1/(4*pi)./sqrt( (-yq(2)-hj*uj(2)/sqrt(uj(2)^2+uj(3)^2)+xj(2)+uj(3).*x.*sin(y)./sqrt(uj(2)^2+uj(3)^2) ).^2+(-yq(1)+hj*uj(1)*uj(3)/sqrt(uj(2)^2+uj(3)^2)+xj(1)+uj(1)*uj(2).*x.*sin(y)./sqrt(uj(2)^2+uj(3)^2)+sqrt(uj(2)^2+uj(3)^2).*x.*cos(y)).^2+(-yq(3)+hj*uj(1)*uj(3)/sqrt(uj(2)^2+uj(3)^2)+xj(1)+uj(1)*uj(2).*x.*sin(y)./sqrt(uj(2)^2+uj(3)^2)+sqrt(uj(2)^2+uj(3)^2).*x.*cos(y)).^2);

if(index==0) %H
    if(n==0)
        f = @(x,y) x.*besselj(0,RB(1,m).*x).*Fai(x,y);
        v = integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(1,RB(1,m))^2)*exp(-RB(1,m)*hj));
    else
        f = @(x,y) x.*besselj(n,RB(n+1,m).*x).*cos(n.*y).*Fai(x,y);
        v = 2*integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(n+1,RB(n+1,m))^2)*exp(-RB(1,m)*hj));
    end
end
if(index==1) %H1
    if(n==0)
        v = 0;
    else
        f = @(x,y) x.*besselj(n,RB(n+1,m).*x).*sin(n.*y).*Fai(x,y);
        v = 2*integral2(f,0,rj,0,2*pi)/(pi*(a^2)*(besselj(n+1,RB(n+1,m))^2)*exp(-RB(1,m)*hj));
    end
end