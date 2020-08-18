function [u, mean_resp] = SQ_exp(rho, d, pi0, N, lam)
%

if ~exist('N','var')
    N=200;
end

if ~exist('lam','var')
    lam=rho;
end

u=zeros(1,N);
u(1)=1; u(2)=rho;
for k=2:N-1
    u(k+1)=rho^((d^k-1)/(d-1))*pi0^((d^(k-1)-1)/(d-1));
end
if lam < 10^(-8)
    mean_resp=1;
else
    mean_resp=(1/lam)*(sum(u(2:end)));
end

end