% real < 0 (take in account nan)
function b=isrealleft(x)
    b=0;
    if(isnan(x))
        b = 0;
        return;
    end
    if(imag(x)!=0)
        b=0;
        return;
    end
    if(real(x)<0)
        b=1;
        return;
    else
        b=0;
        return
    end
    return
end
