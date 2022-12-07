%this function generates a bash script to iterate over the specified
%conditions and call the simulation funtion
clear
datestr=%prefix giving date
SimSetName=%name simulation set 
simfuncname='MergeFunFINAL'; %function named here
outdir=%specify out directory
foldername=[datestr SimSetName];
jobsetnum=1;
replicates=20;
vf=.05; %10.^(-2:.1:-.01); %.05;
alpha=[1];
Dcoefexp=[-4:.5:-2]; %[-4 -3];%[1];%-2:.5:2;
stepsize=[1e-2]; %10.^(-2);%(-3.5:.25:0); %[.01:.01:.1];
Vmin=0.9;
Vmaxstart=1.1;
condlist=combvec(stepsize,alpha,Dcoefexp,vf)';
interval=0.1;
LatticeFlag=0;
%%
%injection no prerun 
prerunflag=0;
Ntot=1500;
Nstart=0;
J=[5 1 1/5 1/10 1/25 1/50 1/100]';
%J=[5 1/100]';
injectionstart=ones(size(J));
injectionend=((Ntot-Nstart)./J);
T=injectionend+2;
%%
%fast nucleation sims with no prerun
prerunflag=0;
Ntot=1000;
Nstart=1000;
T=10^5+2;  
J=0;
injectionstart=T+2;
injectionend=T+2;
%%

hrnum=24.*ones(size(T));
threadnum=20;
Ntot=(injectionend-injectionstart+1).*J+Nstart;
injcondlist=[J T hrnum injectionstart injectionend];
%%
localbatchdest=%insert local batch destination 
mkdir(localbatchdest);
clusterbatchdirec=['/home/dswlee/' foldername '/'];
tmpdirname=['/scratch/gpfs/dswlee/'];
batchsubmitfid=fopen([localbatchdest 'batch_submit.sh'],'w');
fprintf(batchsubmitfid,['mkdir ' outdir '\n']);
for j = 1 :jobsetnum   
for k=1:size(injcondlist,1)
    for i = 1 :size(condlist,1)
    setnum=j;
fnamestr=[SimSetName num2str(i) 'N' num2str(j) '_' num2str(k)];
fid=fopen([localbatchdest fnamestr '.m'],'w');
fprintf(fid, ['clear\noutdir =''' outdir ''';\nreplicates=%d;\n'],replicates);
fprintf(fid,['Nstart=%d;\n'],Nstart);
fprintf(fid,['T=%g;\n'],injcondlist(k,2));
fprintf(fid,['vf=%g;\n'],condlist(i,4));
fprintf(fid,['Vmin=%g;\n'],Vmin);
fprintf(fid,['Vmaxstart=%g;\n'],Vmaxstart);
%fprintf(fid,['Rstart=%d;\n'],Rstart);
fprintf(fid,['alpha=%g;\n'],condlist(i,2));
fprintf(fid,['stepsize=%d;\n'],condlist(i,1));
fprintf(fid,['Dcoefexp=%g;\n'],condlist(i,3));
fprintf(fid,['outTs=ceil(10.^[0:%g:log10(T)]);\n'],interval); %logspace ts
%fprintf(fid,['outTs=[0:%g:T];\n'],interval); %linspace ts
fprintf(fid,['LatticeFlag=%d;\n'],LatticeFlag);
fprintf(fid,['injectionstart=%d;\n'],injcondlist(k,4));
fprintf(fid,['injectionend=%d;\n'],injcondlist(k,5));
fprintf(fid,['J=%d;\n'],injcondlist(k,1));
fprintf(fid,['setnum=%d;\n'],setnum);
fprintf(fid,['prerunflag=%d;\n'],prerunflag);
fprintf(fid,['pc = parcluster(''local'');\n']);
fprintf(fid,['pc.JobStorageLocation = strcat(''' tmpdirname ''''  ', getenv(''SLURM_JOB_ID''));\n']);
fprintf(fid,['poolobj=parpool(pc, %d);\n'],threadnum);
fprintf(fid,['cd ''/home/dswlee'';\n']);
fprintf(fid,['parfor rep=1:replicates\n\trng(''shuffle'')\n\t' simfuncname '(replicates,Nstart,J,injectionstart,injectionend,T,vf,stepsize,alpha,Vmin,Vmaxstart,Dcoefexp,outdir,rep,outTs,LatticeFlag,setnum,prerunflag);\nend\ndelete(poolobj);']);
jobname=[SimSetName num2str(i) 'N' num2str(j) '_' num2str(k) '.sh'];
slurmfid=fopen([localbatchdest  fnamestr '.sh'],'w');
fprintf(slurmfid,'#!/bin/sh\n#SBATCH -N 1\n');
fprintf(slurmfid,'#SBATCH -n 1\n');
fprintf(slurmfid,['#SBATCH -c %d\n'],threadnum);
fprintf(slurmfid,['#SBATCH -J ' jobname(1:end-3) '\n']);
fprintf(slurmfid,['#SBATCH -t ' num2str(injcondlist(k,3)) ':00:00\n']);
fprintf(slurmfid,'#SBATCH --mem-per-cpu=20000\n');
%fprintf(slurmfid,'#SBATCH --mem-per-cpu=10G\n');
fprintf(slurmfid,['#SBATCH --output=' jobname(1:end-3) '\n']);
fprintf(slurmfid,'#SBATCH --mail-user=dswlee@princeton.edu\n');
fprintf(slurmfid,'#SBATCH --mail-type=begin\n');
fprintf(slurmfid,'#SBATCH --mail-type=end\n');
fprintf(slurmfid,[ 'mkdir -p ' tmpdirname  '$SLURM_JOB_ID\n']);
fprintf(slurmfid,'module purge\nmodule load matlab/R2019a\n');
fprintf(slurmfid,['srun matlab -nodesktop -nosplash -r "run(''' clusterbatchdirec fnamestr ''')"\n']);
fprintf(slurmfid, ['rm -rf ' tmpdirname  '$SLURM_JOB_ID\n']);
fprintf(batchsubmitfid, ['sbatch ' clusterbatchdirec fnamestr '.sh\n']);
fclose(fid);
fclose(slurmfid);
end
end
end
fclose(batchsubmitfid);
%%
