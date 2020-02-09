function [roui,faii,si] = GCL_cart2cyl(S,i,x)
u = S.u;
ui = u(i,:);
si = x(1)*ui(1)+x(2)*ui(2)+x(3)*ui(3);
roui = sqrt(x(1)^2+x(2)^2+x(3)^2 - si^2);
if(roui==0)
    faii=0;
elseif(x(2)>=0)
    faii = acos(x(1)/roui);
else 
    faii = acos(-x(1)/roui) + pi;
end    