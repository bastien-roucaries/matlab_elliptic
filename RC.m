% MIT License
%
% Copyright (c) 2018 Laboratoire SATIE
% Copyright (c) 2018 Université de Cergy Pontoise
% Copyright (c) 2018 Bastien Roucariès
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
function f=RC(x,y,errtol=1e-3)
    assert(min(size(x)==size(y)),'x and y should have same size');


    % initialization
    f = 0.*x;

    % vectorize of poor
    for ii=1:numel(x)
        f(ii)=RCscalar(x(ii),y(ii),errtol);
    end
end

% special cases
function [s f]=RCspecial(x,y)
    f=0;
    s=1;

    % nan
    if(isnan(x) || isnan(y))
        f=NaN;
        return;
    end

    % not defined
    if(isrealleft(x))
        f=NaN;
        return;
    end

    if(x==0 && y==0)
        f=NaN;
        return
    end

    if(x==y)
        f=1.0/sqrt(x);
        return;
    end

    % infinity (y is dominant)
    if(isinf(abs(y)))
        f=0;
        return;
    end

    if(isinf(abs(x)))
        % DLMF 19.6.15
        if(y==0)
            f=+inf;
            return;
        else
            f=0;
            return;
        end
    end

    % DLMF §19.6.15
    if(x==0)
        if(real(y)<0 && imag(y) == 0)
            f=0;
            return;
        else
            f=0.5*pi/sqrt(y);
            return;
        end
    end

    % DLMF §19.6.15
    if(y==0)
        if(abs(x)>0)
            f=+inf;
            return;
        else
            f=NaN;
            return;
        end
    end
    s=0;
    f=0;
end

function f=RCscalelow(x,y,errtol,nsl,sqrtnsl)
% scale and loss figures in x (downscale)
% see DLMF 19.6.15
% moreover it will use subnormal
    [x y]=num2cell([x y]./nsl){:};
    [s f]=RCspecial(x,y);
    if(s)
        f=f/sqrtnsl;
        return;
    end
    f=RCgen(x,y,errtol)/sqrtnsl;
    return;
end

function f=RCscaleup(x,y,errtol,nsl,sqrtnsl)
% scale and loss figures in x (downscale)
% see DLMF 19.6.15
% moreover it will use subnormal
    [x y]=num2cell([x y].*nsl){:};
    [s f]=RCspecial(x,y);
    if(s)
        f=f*sqrtnsl;
        return;
    end
    f=RCgen(x,y,errtol)*sqrtnsl;
    return;
end

% scalar version
function f=RCscalar(x,y,errtol)
% Argument limits as set by Carlson (use power of two instead
% of 5.0)
    factorloss = 8.0; % aka 3 bits
    nearestsquareloss = 16.0;
    sqrtnearestsquareloss = 4.0;
    LoLim = factorloss * realmin;
    LoLimS = nearestsquareloss * 2 * realmin;
    UpLim = realmax/factorloss;
    UpLimS = realmax/(2*nearestsquareloss);

    [s f]=RCspecial(x,y);
    if(s)
        return;
    end


    % some special case (huge/small)
    % for Rc y dominate against x
    if(abs(y)>UpLim && abs(x)<LoLimS)
        f=RCscalelow(x,y,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end

    if(abs(y)<LoLim && abs(x)>UpLimS)
        f=RCscaleup(x,y,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end

    % huge
    if(abs(y)>UpLim)
        f=RCscalelow(x,y,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end

    if(abs(x)>UpLim)
        f=RCscalelow(x,y,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return
    end
    % small
    if(abs(y)<LoLim)
        f=RCscaleup(x,y,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
    if(abs(x)<LoLim)
        % see DLMF 19.6.15
        f=RCscaleup(x,y,errtol,nearestsquareloss, ...
                     sqrtnearestsquareloss);
        return;
    end
    f=RCgen(x,y,errtol);
end

% compute by
function f=RCgen(x,y,errtol)
    alamb=0;
    ave=0;
    s=0;
    w=0;
    xt=0;
    yt=0;

    % constant
    C1=0.3;
    C2=(1.0/7.0);
    C3=0.375;
    C4=(9.0/22.0);

    % compute cauchy value
    if (y > 0.0)
        xt=x;
        yt=y;
        w=1.0;
    else
        xt=x-y;
        yt = -y;
        w=sqrt(x)/sqrt(xt);
    end
    contloop = 1;
    while (contloop)
        alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
        xt=0.25*(xt+alamb);
        yt=0.25*(yt+alamb);
        ave=(xt+yt+yt)/3.0;
        s=(yt-ave)/ave;
        contloop = (abs(s) > errtol);
    end
    f=w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
end
