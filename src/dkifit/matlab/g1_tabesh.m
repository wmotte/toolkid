function val=g1_tabesh(vals)

% By Umesh R S on 23 Feb 2010

if vals(2)~=vals(3)
    a=sum(vals)^2/(18*vals(2)*(vals(2)-vals(3))^2);
    b=(vals(3)^2 - 3*vals(2)*vals(3))/sqrt(vals(2)*vals(3));
    
    val=a*(2*vals(2) + b);
else
    val=(vals(1)+2*vals(2))^2/(24*vals(2)^2);
end

return;