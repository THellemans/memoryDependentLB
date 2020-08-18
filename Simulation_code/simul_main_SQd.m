function [mean_resp]=simul_main(d,M,memtype,N,lam,jtype,sp1,sp2,nruns,runlen)
%% This code simulates the memory based SQ(d) model
%M is the memory size
%memtype is the type of memory:
% memtype=1 : No memory
% memtype=2 : Minimal Memory
% memtype=3 : Probe for a bit of memory/a lot of memory when M is large
% memtype=4 : Probe for a bit of memory/a lot of memory when M is large
% memtype=5 : JIQ memory type

mu=zeros(2,1); p=0;
if (jtype == 2)
    SCV=sp1;
    f=sp2;
    mu(1)=(SCV+(4*f-1)+((SCV-1)*(SCV-1+8*f*(1-f)))^(1/2))/(2*1*f*(SCV+1));
    mu(2)=(SCV+(4*(1-f)-1)-((SCV-1)*(SCV-1+8*f*(1-f)))^(1/2))/(2*1*(1-f)*(SCV+1));
    p=1*mu(1)*f;
end
mean_resp=zeros(1,nruns);
idles=zeros(1,N,'logical');
nr_jobs=zeros(1,N);
start_recording=runlen/3;
for nrun=1:nruns
    nrun
    wloads=zeros(1,N);
    ctime=0;
    total_response=0;
    nr_responses=0;
    while ctime < runlen
        iat=exprnd(1/(lam*N));
        ctime=ctime+iat;
        lose_job=find( (wloads-iat)<=0 & nr_jobs > 0);
        extra_loss=zeros(size(lose_job));
        for i=1:length(lose_job)
            extra_loss(i)=wloads(lose_job(i))-iat;
        end
        wloads=max(wloads-iat,0);
        for j=1:length(lose_job)
            i=lose_job(j);
            nr_jobs(i)=nr_jobs(i)-1;
            if jtype==1
                jsize=exprnd(1);
            else
                jsize=exprnd(1/mu(1+(rand>p)));
            end
            if nr_jobs(i) >0
                wloads(i)=max(jsize+extra_loss(j),0);
            else
                if sum(idles) < M && memtype==5
                    idles(i)=true;
                end
            end
        end
        
        if all(~idles)
            qids=ceil(rand(1,d)*N);
            [w,indx]=min(nr_jobs(qids));
            if ctime > start_recording
                to_add=0;
                nr_responses=nr_responses+1;
                for i=1:w
                    if jtype==1
                        extra_job=exprnd(1);
                    else
                        extra_job=exprnd(1/mu(1+(rand>p)));
                    end
                    to_add=to_add+extra_job;
                end
                if w==0
                    to_add=exprnd(1);
                end
                total_response=(nr_responses-1)/(nr_responses)*total_response+(to_add+wloads(qids(indx)))/nr_responses;
            end
            nr_jobs(qids(indx))=w+1;
            if w==0
                if jtype==1
                    jsize=exprnd(1);
                else
                    jsize=exprnd(1/mu(1+(rand>p)));
                end
                wloads(qids(indx))=jsize;
            end
            if (memtype == 2 || memtype == 3 || memtype==4)
                for a=1:length(qids)
                    i=qids(a);
                    if sum(idles)<M && nr_jobs(i)==0
                        idles(i)=true;
                    end
                end
            end
        else
            server=find(idles,1,'first');
            if jtype==1
                jsize=exprnd(1);
            else
                jsize=exprnd(1/mu(1+(rand>p)));
            end
            if ctime > start_recording
                nr_responses=nr_responses+1;
                total_response=(nr_responses-1)/(nr_responses)*total_response+(jsize)/nr_responses;
            end
            wloads(server)=jsize;
            nr_jobs(server)=1;
            idles(server)=false;
            if (memtype==3 || memtype==4)
                qids_extra=ceil(rand(1,d)*N);
                extra_idles=find(wloads(qids_extra)==0);
                mem_size=sum(idles);
                for i=1:length(extra_idles)
                    if mem_size+i<=M
                        idles(qids_extra(extra_idles(i)))=true;
                    else
                        break
                    end
                end
            end
        end
    end
    mean_resp(nrun)=total_response;
end
end