function val=g2_tabesh(vals)

% By Umesh R S on 23 Feb 2010

if vals(2)~=vals(3)
    a=sum(vals)^2/(3*(vals(2)-vals(3))^2);
    b=(vals(2)+vals(3))/sqrt(vals(2)*vals(3));
    
    val=a*(b -2);
else
    val=6*(vals(1)+2*vals(2))^2/(72*vals(2)^2);
end

return;