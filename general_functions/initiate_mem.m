function [pi0,probesused] = initiate_mem( memtype, rho, d, A )
% memtype = 1 => No Memory
% memtype = 2 => IP
% memtype = 3 => CP
% memtype = 4 => BCP
% memtype = 5 => ISM

if rho == 0 && memtype ~= 1
    pi0=0;
    probesused=0;
    return
end

if memtype==1 %no memory
    pi0=1;
elseif memtype==2 %IP
    pi0=1/(rho^d+(1-rho)*d);
elseif memtype==3 %CP
    if d<1/(1-rho)
        pi0=(1-d*(1-rho))/rho^d;
    else
        pi0=0;
    end
elseif memtype==4 %BCP
    if d==2
        if abs(rho-1/2)<10^(-4)
            pi0=1/(1+A);
        else
            pi0=(1-((1-rho)/rho)^2) / (1-((1-rho)/rho)^(2*(A+1)));
        end
    else
        [M] = define_M_BCP(rho, d, A);
        pis=M^1000;
        pi0=pis(1);
    end
elseif memtype==5 %ISM
    pi0=(1-(1-rho^d)^(1/(A+1)))/(rho^d);
    probesused=(1-pi0)+pi0*(1-rho^d)/(1-rho)+pi0*(1-rho^d);
end
if  memtype~=5
    probesused=(1-pi0*rho^d)/(1-rho);
end
if memtype==3
    probesused=d;
end

end