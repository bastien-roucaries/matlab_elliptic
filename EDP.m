% D(phi,k)
function f=EDP(phi,ak)
    cc=(cos(phi))^2;
    s=sin(phi);
    q=(1.0-s*ak)*(1.0+s*ak);
    f=s^3*RD(cc,q,1.0)/3.0;
end