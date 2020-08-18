addpath('../general_functions')
addpath('../../Chapter 1')

%% Figure 1 : SQ(5), memory size = 4 ISM memory & exp(1) job sizes.
d=5;
A=4;
lam_iter=linspace(0.01,1,10^3);
memtypes={"No Memory", "IP", "CP", "BCP", "ISM"};
ER=zeros(length(memtypes), length(lam_iter));
probesused=zeros(size(ER));
pi0_holder=zeros(size(ER));
N=200;

for i=1:length(memtypes)
    for j=1:length(lam_iter)
        lam=lam_iter(j);
        [pi0_holder(i,j),probesused(i,j)] = initiate_mem( i, lam, d, A );
        [~, ER(i,j)] = SQ_exp(lam, d, pi0_holder(i,j), N, lam);
    end
end

save('../data/figure1.mat', 'pi0_holder', 'probesused', 'ER', 'lam_iter', 'd', 'A', 'memtypes')
