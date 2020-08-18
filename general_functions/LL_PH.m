function [Fbar, FRbar, ER] = LL_PH(rho, d, lam, pi0, Gbar, g, alpha, A, w_range)
%

if length(alpha) == 1
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
            if temp<10^(-8)
                break
            end
            n=n+1;
        end
    end
else
    mu=-A*ones(length(A),1);
    Fbar=zeros(size(w_range));
    Fbar(1)=rho;
    dw=w_range(2)-w_range(1);
    xi=zeros(length(A), length(w_range));
    for n=1:length(w_range)-1
        Fbar(n+1)=Fbar(n)-lam*dw*(Gbar(n)+pi0*(-Fbar(n)^d+alpha*xi(:,n)));
        xi(:,n+1)=xi(:,n)+dw*(A*xi(:,n)+Fbar(n)^d*mu);
    end
    FRbar=(1-pi0)*Gbar+pi0*(Gbar+convolution(Fbar.^d, g, dw));
    ER=trapz(FRbar)*dw;
end
end