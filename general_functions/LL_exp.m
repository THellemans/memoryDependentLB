function [Fbar, FRbar, ER] = LL_exp(rho, d, pi0, w_range, tol)
%


Fbar=(rho*pi0+(rho^(1-d)-rho*pi0)*exp((d-1)*w_range)).^(1/(1-d));
FRbar=(rho^d*pi0+(1-rho^d*pi0)*exp((d-1)*w_range)).^(1/(1-d));
if d==2
    ER=-log(1-pi0*rho^2)/(pi0*rho^2);
else
    ER=0;
    n=0;
    while true
        temp=pi0^n*rho^(d*n)/(1+n*(d-1));
        ER=ER+temp;
        if temp<tol
            break
        end
        n=n+1;
    end
end

end