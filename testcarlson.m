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
function testcarlson
    % RC
    tol(RC(0,1.0/4.0),pi)
    tol(RC(9.0/4.0,2.0),log(2))
    tol(RC(0,i),(1-i)*1.1107207345396)
    tol(RC(-i,i), 1.2260849569072-i*0.34471136988768)
    tol(RC(1.0/4.0,-2.0),log(2)/3.0)
    tol(RC(i,-1.0),0.77778596920447+i*0.19832484993429)
    assert(isnan(RC(NaN,1.0)))
    assert(isnan(RC(1.0,NaN)))
    assert(isnan(RF(NaN,0,0)))
    % RF
    tol(RF(1,2,0),1.3110287771461)
    tol(RF(i,-i,0),1.8540746773014)
    tol(RF(0.5,1.0,0),1.8540746773014)
    tol(RF(2,3,4),0.58408284167715)
    tol(RF(i,-i,2),1.0441445654064)
    tol(RF(i-1,i,1-i), 0.93912050218619-i*0.53296252018635)
    % RD
    tol(RD(0.0,2.0,1.0), 1.7972103521034)
    tol(RD(2.0,3.0,4.0),0.16510527294261)
    tol(RD(i,-i,2.0),0.65933854154220)
    tol(RD(0,i,-i),1.2708196271910+i*2.7811120159521)
    tol(RD(0,i-1,i),-1.8577235439239- i*0.96193450888839)
    tol(RD(-2-i,-i,-1+i),1.8249027393704 -i*1.2218475784827)
    % RJ
    tol(RJ(0,1,2,3),0.77688623778582)
    tol(RJ(2,3,4,5),0.14297579667157)
    % Does not work
    %tol(RJ(2,3,4,-1+i),0.13613945827771-i*0.38207561624427)%
    tol(RJ(i,-i,0,2), 1.6490011662711)
    tol(RJ(-1+i,-1-i,1,2),0.94148358841220)
    tol(RJ(i,-i, 0,1-i),1.8260115229009+i*1.2290661908643)
    % does not work
    %tol(RJ(-1 + i, -1+i,1,-3+i),-0.61127970812028 -i*1.0684038390007)%
    tol(RJ(-1+i,-2-i,-i,-1+i),1.8249027393704-i*1.2218475784827)
    tol(RJ(2, 3,4,-0.5),0.24723819703052)
    tol(RJ(2, 3,4, -5),-0.12711230042964)
    % from AS p608
    tol(EF(pi/2,0),pi/2)
    tol(EF(pi/2,sqrt(0.2)),1.659623598610528)
    tol(EF(pi/2,sqrt(0.48)),1.837491363355796)
    % from AS p609
    %tol(EE(sqrt(0.99)),1.015993546)
    % D(k)
    tol(EDP(pi/3,0.5),(EF(pi/3,0.5)-EEP(pi/3,0.5))/0.5^2)
    tol(ED(0.5),(EK(0.5)-EE(0.5))/0.5^2)
    % from mathematica
    % from mathematica
    tol(EPi(0.4,sqrt(0.6)), ...
        2.590921156555220293067792172295651236897020458546480266071)
     tol(EPiP(pi/2,0.4,sqrt(0.6)), ...
        2.590921156555220293067792172295651236897020458546480266071)
    % EllipticPi[0.4,pi/3,0.6] !!! see mathematica convention
    tol(EPiP(pi/3,0.4,sqrt(0.6)), ...
             1.353572359896644385640544495988724325661999937255699041262)
    
end

function tol(x,c,t=eps*1000)
   test=abs(x-c)/abs(c);
   assert(test<t,sprintf('value is outside tol: get %e tol %e',test,t));
end

