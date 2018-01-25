% F(phi,k)
function f=EF(phi,ak)
   s=sin(phi);
   f=s.*RF(cos(phi).^2,(1.0-s.*ak).*(1.0+s.*ak),1.0);
end