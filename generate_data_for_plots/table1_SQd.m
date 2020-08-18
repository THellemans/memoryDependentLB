addpath('../general_functions')
addpath('../../Chapter 1')
addpath('../Simulation_code')

qbound=1000;
runtime_times_N=10^6;
N_iter=[10 20 50 100 200 500 1000 3000];
%% d=4, memtype 2, lam=0.9, jtype=exp.
M=10000; lam=0.9; jtype=1; sp1=0; sp2=0; nruns=40;
d=4; memtype=2;

for N=N_iter
    runlen=runtime_times_N/N;
    [mean_resp]=simul_main_SQd_mex(d,M,memtype,N,lam,jtype,sp1,sp2,nruns,runlen);
    save(char(strcat("../data/data_SQd_memtype2_N",int2str(N))), 'lam','memtype','M','N','mean_resp');
end
[pi0,~] = initiate_mem( memtype, lam, d, M );
[~, ER] = SQ_exp(lam, d, pi0, qbound, lam);
save(char(strcat("../data/data_SQd_memtype2_Ninfty.mat")), 'lam','memtype','d','ER','pi0');

%% d=3, memtype 3, lam=0.8, jtype=exp.
M=10000; lam=0.8; jtype=1; sp1=0; sp2=0; nruns=40;
d=3; memtype=3;

for N=N_iter
    runlen=runtime_times_N/N;
    [mean_resp]=simul_main_SQd_mex(d,M,memtype,N,lam,jtype,sp1,sp2,nruns,runlen);
    save(char(strcat("../data/data_SQd_memtype3_N",int2str(N))), 'lam','memtype','M','N','mean_resp');
end
[pi0,~] = initiate_mem( memtype, lam, d, M );
[~, ER] = SQ_exp(lam, d, pi0, qbound, lam);
save(char(strcat("../data/data_SQd_memtype3_Ninfty.mat")), 'lam','memtype','d','ER','pi0');

%% d=3, memtype 4, lam=0.8, M=5, jtype=HExp.
M=5; lam=0.8; jtype=2; sp1=2; sp2=1/2; nruns=40;
d=3; memtype=4;

for N=N_iter
    runlen=runtime_times_N/N;
    [mean_resp]=simul_main_SQd_mex(d,M,memtype,N,lam,jtype,sp1,sp2,nruns,runlen);
    save(char(strcat("../data/data_SQd_memtype4_N",int2str(N))), 'lam','memtype','M','N','mean_resp');
end

[pi0,~] = initiate_mem( memtype, lam, d, M );
[~, ER] = SQ_PH(lam, d, 200, sp1, pi0);
save(char(strcat("../data/data_SQd_memtype4_Ninfty.mat")), 'lam','memtype','d','ER','pi0');

%% d=2, memtype 5, lam=0.85, M=10, jtype=HExp.
M=10; lam=0.85; jtype=2; sp1=3; sp2=1/2; nruns=40;
d=2; memtype=5;

for N=N_iter
    runlen=runtime_times_N/N;
    [mean_resp]=simul_main_SQd_mex(d,M,memtype,N,lam,jtype,sp1,sp2,nruns,runlen);
    save(char(strcat("../data/data_SQd_memtype5_N",int2str(N))), 'lam','memtype','M','N','mean_resp');
end

[pi0,~] = initiate_mem( memtype, lam, d, M );
[~, ER] = SQ_PH(lam, d, 200, sp1, pi0);
save(char(strcat("../data/data_SQd_memtype5_Ninfty.mat")), 'lam','memtype','d','ER','pi0');