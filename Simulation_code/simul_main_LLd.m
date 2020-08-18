function [mwloads, mean_resp]=simul_main(d,M,memtype,N,lam,jtype,sp1,sp2,nruns,runlen)
% M is the memory size
% memtype is the type of memory:
% memtype=1 : No memory
% memtype=2 : Minimal Memory
% memtype=3 : Probe for a bit of memory/a lot of memory when M is large
% memtype=4 : (here same as 3)
% memtype=5 : JIQ like memory

mu=zeros(2,1); p=0;
if (jtype == 2)
    SCV=sp1;
    f=sp2;
    mu(1)=(SCV+(4*f-1)+((SCV-1)*(SCV-1+8*f*(1-f)))^(1/2))/(2*1*f*(SCV+1));
    mu(2)=(SCV+(4*(1-f)-1)-((SCV-1)*(SCV-1+8*f*(1-f)))^(1/2))/(2*1*(1-f)*(SCV+1));
    p=1*mu(1)*f;
end
mwloads=zeros(1,nruns);
idles=zeros(1,N,'logical');
mean_resp=zeros(1,nruns);
start_recording=runlen/3;
for nrun=1:nruns
	nrun
	total_response=0;
    nr_responses=0;
    wloads=zeros(1,N);
    ctime=0; octime=0;
    nm=0; mwload=0;
    while ctime < runlen
        iat=exprnd(1/(lam*N));
        ctime=ctime+iat;
        if memtype==5
            become_idle=find( ((wloads-iat)<=0) & (wloads>0) );
            for i=1:length(become_idle)
                if sum(idles) == M
                    break
                end
                idles(become_idle(i))=true;
            end
        end
        wloads=max(wloads-iat,0);
        
        if jtype==1
            jsize=exprnd(1);
        else
            jsize=exprnd(1/mu(1+(rand>p)));
        end
        
        if all(~idles)
            qids=ceil(rand(1,d)*N);
            [w,indx]=min(wloads(qids));
            wloads(qids(indx))=w+jsize;
            if ctime > start_recording
                nr_responses=nr_responses+1;
                total_response=(nr_responses-1)/(nr_responses)*total_response+(w+jsize)/nr_responses;
            end
            if (memtype == 2 || memtype == 3 || memtype == 4)
                idles(qids(wloads(qids)==0))=true;
            end
        else
        	if ctime > start_recording
                nr_responses=nr_responses+1;
                total_response=(nr_responses-1)/(nr_responses)*total_response+(jsize)/nr_responses;
            end
            server=find(idles,1,'first');
            wloads(server)=jsize;
            idles(server)=false;
            if (memtype==3 || memtype == 4)
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
        if (ctime > runlen/3 && ctime-octime > N*runlen/200000)
            mwload=mwload+sum(wloads)/N;
            octime=ctime;
            nm=nm+1;
        end
    end
    mean_resp(nrun)=total_response;
    mwloads(nrun)=mwload/nm;
end
end