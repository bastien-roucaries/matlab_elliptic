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