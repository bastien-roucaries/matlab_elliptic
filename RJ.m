function f=RJ(x, y, z, p, errtol=1e-4)
    assert(min([min(size(x)==size(y)),min(size(x)==size(z)),min(size(x)==size(p))]), ...
       'x y z size should be equal');

    f=0.*x;

    for ii=1:numel(x)
        if(isnan(x(ii)) || isnan(y(ii)) || isnan(z(ii)) || isnan(p(ii)))
            f(ii)=NaN;
            continue
        end
        % RJ is symetric sort by magnetude
        s=sort([x(ii),y(ii) z(ii)]);
        f(ii)=RJscalar(s(1),s(2),s(3),p,errtol);
    end
end

function [s f]=RJspecial(x,y,z,p)
    s=1;
    if(isnan(x) || isnan(y) || isnan(z) || isnan(p))
        f=NaN;
        return
    end
    % RJ is symetric sort by magnetude
    [xv]=sort([x y z]);
    x=xv(1);
    y=xv(2);
    z=xv(3);
    if(isrealleft(x) || isrealleft(y) ||  isrealleft(z))
        f=NaN;
        return
    end
    
    % 0
    if(x==0 && y ==0)
        f=inf;
        return
    end
    
    % infinity
    if(isinf(abs(x)) || isinf(abs(y)) || isinf(abs(z)) || isinf(abs(p)))
        f=0;
    end
    
    % equal
    if(x==y && x==z && z==p)
        if(x*sqrt(x)==0)
            f=+inf;
        else
            f=1.0/(x*sqrt(x));
        end
        return;
    end
    
    %RD
    if(z==p)
        f=RD(x,y,z);
        return;
    end
    s=0;
    f=0;
 end

function f=RJscalelow(x,y,z,p,errtol,nsl,sqrtnsl)
% scale and loss figures in x (downscale)
% see DLMF 19.6.15
% moreover it will use subnormal
    [x y z p]=num2cell([x y z p]./nsl){:};
    [s f]=RDspecial(x,y,z,p);
    if(s)
        f=f/sqrtnsl;
        return;
    end
    f=RJgen(x,y,z,p,errtol)/sqrtnsl;
    return;    
end

function f=RJscaleup(x,y,z,p,errtol,nsl,sqrtnsl)
% scale and loss figures in x (downscale)
% see DLMF 19.6.15
% moreover it will use subnormal
    [x y z p]=num2cell([x y z p].*nsl){:};
    [s f]=RJspecial(x,y,z,p);
    if(s)
        f=f*sqrtnsl;
        return;
    end
    f=RJgen(x,y,z,p,errtol)*sqrtnsl;
    return;    
end

 

function f=RJscalar(x,y,z,p,errtol)
   [s f]=RJspecial(x,y,z,p);
   if(s)
       return
   end
   
   % Argument limits as set by Carlson (use power of two instead
   % of 5.0)
   factorloss = 16.0; % aka 4 bits
   nearestsquareloss = 16.0; 
   sqrtnearestsquareloss = 4.0*16;
   LoLim = factorloss * realmin;
   LoLimS = nearestsquareloss * 2 * realmin;
   UpLim = realmax/factorloss;
   UpLimS = realmax/(2*nearestsquareloss);
   
   % some special case (huge/small)
    % for Rc y dominate against x
    if(abs(p)>UpLim && abs(x)<LoLimS && abs(y)<LoLimS && abs(p)<LoLimS)
        f=RJscalelow(x,y,z,p,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    
    if(abs(p)<LoLim && abs(x)>UpLimS && abs(y)>UpLimS && abs(p)>UpLimS)
        f=RJscaleup(x,y,z,p,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
            
    % huge
    if(abs(p)>UpLim)
        f=RJscalelow(x,y,z,p,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    
    if(abs(x)>UpLim || abs(y) > UpLim || abs(z) > UpLim)
        f=RJscalelow(x,y,z,p,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    % small
    if(abs(p)<LoLim)
        f=RJscaleup(x,y,z,p,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
    if(abs(x)<LoLim || abs(y)<LoLim || abs(z) < LoLim)
        % see DLMF 19.6.15
        f=RJscaleup(x,y,z,p,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
   
   f=RJgen(x,y,z,p,errtol);
end



    
function f=RJgen(x,y,z,p,errtol)
    C1=(3.0/14.0);
    C2=(1.0/3.0);
    C3=(3.0/22.0);
    C4=(3.0/26.0);
    C5=(0.75*C3);
    C6=(1.5*C4);
    C7=(0.5*C2);
    C8=(C3+C3);
    sum=0.0;
    fac=1.0;
    if (p > 0.0)
        xt=x;
        yt=y;
        zt=z;
        pt=p;
    else
        s=sort([x y z]);
        xt=s(1);
        yt=s(2);
        zt=s(3);
        a=1.0/(yt-p);
        b=a*(zt-yt)*(yt-xt);
        pt=yt+b;
        rho=xt*zt/yt;
        tau=p*pt/yt;
        rcx=RC(rho,tau);
    end
    cont=1;
    while(cont)
        sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)^2;
        beta=pt*(pt+alamb)^2;
        sum += fac*RC(alpha,beta);
        fac=0.25*fac;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        pt=0.25*(pt+alamb);
        ave=0.2*(xt+yt+zt+pt+pt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
        delp=(ave-pt)/ave;
        cont=(max(abs([delx dely delz delp])) > errtol);
    end
    ea=delx*(dely+delz)+dely*delz;
    eb=delx*dely*delz;
    ec=delp*delp;
    ed=ea-3.0*ec;
    ee=eb+2.0*delp*(ea-ec);
    ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+...
                     delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
    if (p <= 0.0) 
        ans=a*(b*ans+3.0*(rcx-RF(xt,yt,zt)));
    end
    f=ans;
end
