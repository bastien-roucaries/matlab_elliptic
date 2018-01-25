function f=RF(x, y, z, errtol=1e-3)
    assert(min([min(size(x)==size(y)),min(size(x)==size(z))]), ...
       'x y z size should be equal');

    f=0.*x;

    for ii=1:numel(x)
        if(isnan(x(ii)) || isnan(y(ii)) || isnan(z(ii)))
            f(ii)=NaN;
            continue
        end
        % Rf is symetric sort by magnetude
        s=sort([x(ii),y(ii),z(ii)]);
        f(ii)=RFscalar(s(1),s(2),s(3),errtol);
    end
end

function [s f]=RFspecial(x,y,z)
    s=1;
    if(isnan(x) || isnan(y) || isnan(z))
        f=NaN;
        return
    end
    % Rf is symetric sort by magnetude
    [xv]=sort([x y z]);
    x=xv(1);
    y=xv(2);
    z=xv(3);
    if(isrealleft(x) || isrealleft(y) ||  isrealleft(z))
        f=NaN;
        return
    end
    
    % 0
    if(x==0 && y ==0 && z==0)
        f=NaN;
        return
    end
    
    % equal
    if(x==y && x==z)
        f=1.0/sqrt(x);
        return;
    end
    
    % RC
    if(x==y)
        f=RC(z,y);
        return
    end
    if(x==z)
        f=RC(y,z);
        return
    end
    if(y==z)
        f=RC(x,y);
        return
    end
    
    % infinity
    if(isinf(abs(x)) || isinf(abs(y)) || isinf(abs(z)))
        f=0;
    end
    s=0;
    f=0;
    end


function f=RFscalar(x,y,z,errtol)
   [s f]=RFspecial(x,y,z);
   if(s)
       return
   end
   
   % Argument limits as set by Carlson (use power of two instead
   % of 5.0)
   factorloss = 16.0; % aka 4 bits
   nearestsquareloss = 16.0; 
   sqrtnearestsquareloss = 4.0;
   LoLim = factorloss * realmin;
   LoLimS = nearestsquareloss * 2 * realmin;
   UpLim = realmax/factorloss;
   UpLimS = realmax/(2*nearestsquareloss);
   
   % prefer to decrease in order to use subnormal
   if(abs(z) > UpLim)
       % use DLMF 19.20.1
       [x y z]=num2cell([x y z]./nearestsquareloss){:};
       [s f]=RFspecial(x,y,z);
       if(s)
           f=f./sqrtnearestsquareloss;
           return;
       end
       f=RFgen(x,y,z,errtol)/sqrtnearestsquareloss;
       return;
   end
   
   if(abs(x) < LoLim)
       [x y z]=num2cell([x y z].*nearestsquareloss){:};
       [s f]=RFspecial(x,y,z);
       if(s)
           f=f*sqrtnearestsquareloss;
           return;
       end
       % use DLMF 19.20.1
       f=sqrtnearestsquareloss*RFgen(x,y,z,errtol);
       return;
   end 
   f=RFgen(x,y,z,errtol);
end



    
function f=RFgen(x,y,z,errtol)
    sqrtx = 0;
    sqrty = 0;
    sqrtz = 0;
    alamb = 0;
    ave = 0;
    delx = 0;
    dely= 0;
    delz = 0;
    epslon = 0;
    e2 = 0;
    e3 = 0;
 
    % constant
    C1 = 1.0/24.0;
    C2 = 0.1;
    C3 = 3.0/44.0;
    C4 = 1.0/14.0;

    % general case see numerical reciptes in C p264
    xt = x;
    yt = y;
    zt = z;
    contloop = 1;
    while (contloop)
        sqrtx=sqrt(xt);
        sqrty=sqrt(yt);
        sqrtz=sqrt(zt);
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        zt=0.25*(zt+alamb);
        ave=(xt+yt+zt)/3.0;
        delx=(ave-xt)/ave;
        dely=(ave-yt)/ave;
        delz=(ave-zt)/ave;
        epslon = max([abs(delx),abs(dely),abs(delz)]);
        contloop = epslon >= errtol;
    end
    e2=delx*dely-delz*delz;
    e3=delx*dely*delz;
    f=(1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
end
