function val=f1_tabesh(vals)

% By Umesh R S on 23 Feb 2010

errtol=1e-6;

if (vals(1)~=vals(2) && vals(1)~=vals(3))
    a=sum(vals)^2/(18*(vals(1)-vals(2))*(vals(1)-vals(3)));
    b=(3*vals(1)^2 - vals(1)*vals(2) - vals(1)*vals(3) - vals(2)*vals(3))/(3*vals(1)*sqrt(vals(2)*vals(3)));
    c=sqrt(vals(2)*vals(3))/vals(1);
    
    val=a*(c*rf(vals(1)/vals(2),vals(1)/vals(3),1,errtol) + b*rd(vals(1)/vals(2),vals(1)/vals(3),1,errtol) -1);
    
elseif vals(1)==vals(2) && vals(2)==vals(3)
    val=0.2;
    
elseif vals(1)==vals(2)
    val=0.5*f2_tabesh([vals(3) vals(1) vals(1)]);
    
elseif vals(1)==vals(3)
    val=0.5*f2_tabesh([vals(2) vals(1) vals(1)]);
    
end

return;