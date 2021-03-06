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
% see B. C. Carlson (1995) Numerical computation of real or complex
% elliptic integrals. Numer. Algorithms 10 (1-2), pp. 13–26.
% DLMF 19.25.14 	
function f=EPiP(phi,a2,ak)
    cc=(cos(phi))^2;
    s=sin(phi);
    q=(1.0-s.*ak).*(1.0+s.*ak);
    f=s^3*a2.*RJ(cc,q,1.0,1.0-a2.*s.^2)/3.0+EF(phi,ak);
end