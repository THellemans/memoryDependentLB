addpath('../general_functions')

load('../data/figure1.mat')

%% figure 1a
x=lam_iter; y=ER; legend_entries=memtypes;
legend_loc='northwest';
xLab='$\lambda$'; yLab='$E[R]$'; xLim=[0.5, 1]; yLim=[0, 5]; filename="../figures/fig1a.pdf"; logplot=false;
make_and_save_fig(x, y, legend_entries, legend_loc, xLab, yLab, xLim, yLim, filename, logplot);

%% figure 1b
x=lam_iter; y=probesused; legend_entries=memtypes;
legend_loc='west';
xLab='$\lambda$'; yLab='Probes Used'; xLim=[0, 1]; yLim=[0, 5.2]; filename="../figures/fig1b.pdf"; logplot=false;
make_and_save_fig(x, y, legend_entries, legend_loc, xLab, yLab, xLim, yLim, filename, logplot);

%% figure 1c
x=lam_iter; y=pi0_holder; legend_entries=memtypes;
legend_loc='northwest';
xLab='$\lambda$'; yLab='$\pi_0$'; xLim=[0, 1]; yLim=[0, 1]; filename="../figures/fig1c.pdf"; logplot=false;
make_and_save_fig(x, y, legend_entries, legend_loc, xLab, yLab, xLim, yLim, filename, logplot);
