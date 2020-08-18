function [M] = define_M_BCP(rho, d, A)
%

M=zeros(A+1,A+1);
M(1,1)=rho^d+nchoosek(d,1)*rho^(d-1)*(1-rho);
for ell=2:d
    M(1,ell)=nchoosek(d,ell)*rho^(d-ell)*(1-rho)^ell;
end
for k=2:A+1
    for ell=(k-1):min((d+k-1),A+1)
        M(k,ell)=nchoosek(d,ell+1-k)*rho^(d-(ell+1-k))*(1-rho)^(ell+1-k);
    end
end
for k=(A-d+3):(A+1)
    M(k,A+1)=0;
    for j=(A-k+2):d
        M(k,A+1)=M(k,A+1)+nchoosek(d,j)*rho^(d-j)*(1-rho)^j;
    end
end
end