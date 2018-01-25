% E(phi,k)
function f=EEP(phi,ak)
    assert(min(size(phi)==size(ak)), ...
       'phi and ak size should be equal');
    
    f=0.*phi;
% see remark in DLMF §19.36(i)
    for ii=1:numel(ak)
        f(ii)=EEPscalar(phi(ii),ak(ii));
    end
end

function f=EEPscalar(phi,ak)
    % see remark in DLMF §19.36(i)
    if(real(ak^2)>0.5 && real(phi) > pi/4 && 0)
        s=sin(phi);
        cc=(cos(phi))^2;
        kp2=1-ak^2;
        q=(1.0-s*ak)*(1.0+s*ak)
        RD(pi/2,0,q)
        f=s*(kp2*RF(cc,q,1.0)+((s*ak)^2)*RD(cc,1.0,q)/3.0)+...
          ak*sqrt((cc*s^2)/(1-s^2*ak^2))
    else
        s=sin(phi);
        cc=(cos(phi))^2;
        q=(1.0-s*ak)*(1.0+s*ak);
        f=s*(RF(cc,q,1.0)-((s*ak)^2)*RD(cc,q,1.0)/3.0);
    end
end