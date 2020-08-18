addpath('../data')

ER_holder=zeros(8,9);
N_vals={'10','20','50','100','200','500','1000','3000','infty'};
for memtype=2:5
    for j=1:length(N_vals)
        N_str=N_vals{j};
        full_filename=strcat('data_LLd_memtype',int2str(memtype),'_N',N_str,'.mat');
        load(strcat('data/',full_filename));
        if strcmp(N_str,'infty')
            ER_holder(memtype-1,j)=ER;
        else
            ER_holder(memtype-1,j)=mean(mean_resp);
        end
        full_filename=strcat('data_SQd_memtype',int2str(memtype),'_N',N_str,'.mat');
        load(strcat('data/',full_filename));
        if strcmp(N_str,'infty')
            ER_holder(memtype+3,j)=ER;
        else
            ER_holder(memtype+3,j)=mean(mean_resp);
        end
    end
end

ER_holder