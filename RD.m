function f=RD(x, y, z, errtol=1e-4)
    assert(min([min(size(x)==size(y)),min(size(x)==size(z))]), ...
       'x y z size should be equal');

    f=0.*x;

    for ii=1:numel(x)
        if(isnan(x(ii)) || isnan(y(ii)) || isnan(z(ii)))
            f(ii)=NaN;
            continue
        end
        % RD is symetric sort by magnetude
        s=sort([x(ii),y(ii)]);
        f(ii)=RDscalar(s(1),s(2),z,errtol);
    end
end

function [s f]=RDspecial(x,y,z)
    s=1;
    if(isnan(x) || isnan(y) || isnan(z))
        f=NaN;
        return
    end
    % Rd is symetric sort by magnetude
    [xv]=sort([x y]);
    x=xv(1);
    y=xv(2);
    if(isrealleft(x) || isrealleft(y) ||  z==0)
        f=NaN;
        return
    end
    
    % 0
    if(x==0 && y ==0)
        if(z==0)
            f=NaN;
        else
            f=+inf;
        end
        return
    end
    
    if(x==0 && y==z)
        if(z*sqrt(z)==0)
            f=inf;
        else
            % DLMF 19.20.1
            f=0.75*pi/(y*sqrt(y));
        end
        return;
    end
    
    % infinity
    if(isinf(abs(x)) || isinf(abs(y)) || isinf(abs(z)))
        f=0;
    end
    
    % equal
    if(x==y && x==z)
        f=1.0/(x*sqrt(x));
        return;
    end
    
    s=0;
    f=0;
 end

function f=RDscalelow(x,y,z,errtol,nsl,sqrtnsl)
% scale and loss figures in x (downscale)
% see DLMF 19.6.15
% moreover it will use subnormal
    [x y z]=num2cell([x y z]./nsl){:};
    [s f]=RDspecial(x,y,z);
    if(s)
        f=f/sqrtnsl;
        return;
    end
    f=RDgen(x,y,z,errtol)/sqrtnsl;
    return;    
end

function f=RDscaleup(x,y,z,errtol,nsl,sqrtnsl)
% scale and loss figures in x (downscale)
% see DLMF 19.6.15
% moreover it will use subnormal
    [x y z]=num2cell([x y z].*nsl){:};
    [s f]=RDspecial(x,y,z);
    if(s)
        f=f*sqrtnsl;
        return;
    end
    f=RDgen(x,y,z,errtol)*sqrtnsl;
    return;    
end

 

function f=RDscalar(x,y,z,errtol)
   [s f]=RDspecial(x,y,z);
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
    if(abs(z)>UpLim && abs(x)<LoLimS && abs(y)<LoLimS)
        f=RDscalelow(x,y,z,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    
    if(abs(z)<LoLim && abs(x)>UpLimS && abs(y)>UpLimS)
        f=RDscaleup(x,y,z,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
            
    % huge
    if(abs(z)>UpLim)
        f=RDscalelow(x,y,z,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    
    if(abs(x)>UpLim || abs(y) > UpLim)
        f=RDscalelow(x,y,z,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    % small
    if(abs(z)<LoLim)
        f=RDscaleup(x,y,z,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
    if(abs(x)<LoLim || abs(z)<LoLim)
        % see DLMF 19.6.15
        f=RDscaleup(x,y,z,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
   
   f=RDgen(x,y,z,errtol);
end



    
function f=RDgen(x,y,z,errtol)
    C1=(3.0/14.0);
    C2=(1.0/6.0);
    C3=(9.0/22.0);
    C4=(3.0/26.0);
    C5=(0.25*C3);
    C6=(1.5*C4);
    
    alamb=0;
    ave=0;
    delx=dely=delz=0;
    ea=eb=ec=ed=ee=0;
    fac=0;
    sqrtx=sqrty=sqrtz=0;
    sum=0;
    xt=yt=zt=0;

    xt=x;
    yt=y;
    zt=z;
    sum=0.0;
    fac=1.0;
    contloop=1;
    while(contloop)
        sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        sum += fac/(sqrtz*(zt+alamb));
        fac=0.25*fac;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        ave=0.2*(xt+yt+3.0*zt);
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
        contloop = max(abs([delx dely delz]) > errtol);
    end
    ea=delx*dely;
    eb=delz*delz;
    ec=ea-eb;
    ed=ea-6.0*eb;
    ee=ed+ec+ec;
    f=3.0*sum+...
      fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
end
