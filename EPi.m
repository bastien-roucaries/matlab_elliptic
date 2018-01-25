% D(phi,k)
function f=EPi(a2,ak)
   f=(a2/3).*RJ(0,1-ak^2,1,1-a2)+EK(ak);
end