function val=f2_tabesh(vals)

% By Umesh R S on 23 Feb 2010
errtol=1e-6;

if vals(2)~=vals(3)
    a=sum(vals)^2/(3*((vals(2)-vals(3))^2));
    b=(2*vals(1)-(vals(2)+vals(3)))/(3*sqrt(vals(2)*vals(3)));
    c=(vals(2)+vals(3))/(sqrt(vals(2)*vals(3)));
    
    val=a*(c*rf(vals(1)/vals(2),vals(1)/vals(3),1,errtol) + b*rd(vals(1)/vals(2),vals(1)/vals(3),1,errtol) -2);
    
elseif vals(3)~=vals(1)
    a=(vals(1)+2*vals(3))^2/(144*vals(3)^2*(vals(1)-vals(3))^2);
    b=vals(3)*(vals(1)+2*vals(3));
    c=(1-(vals(1)/vals(3)));
    d=vals(1)*(vals(1)-4*vals(3));
    
    if c > 0
        c=sqrt(c);
        d=d*(1/c)*atan(c);
    else
        c=sqrt(-c);
        d=d*(1/c)*atan(c);
    end
    
    val=6*a*(b+d);
    
else
    val=6/15;
end

return;