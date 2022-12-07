%this script parses all simulation results in a particular directory and
%saves them to a single matlab workspace which can be passed to
%"workspace_parse_function" to finally read out desired variables

clear
directory=%folder containing csv outputs of simulation

dirobj=dir([directory '*R.csv']);
[fnames{1:length(dirobj)}]=dirobj(:).name;

alphaexp=zeros(length(dirobj),1);
N=zeros(length(dirobj),1);
vf=zeros(length(dirobj),1);
t=zeros(length(dirobj),1);
rep=zeros(length(dirobj),1);

for i = 1:length(fnames)
      filename=fnames{i};
    underlocs=strfind(filename,'_');
    alphaexp(i,1)=str2num(filename(underlocs(1)+1:underlocs(2)-1));
    N(i,1)=str2num(filename(underlocs(3)+1:underlocs(4)-1));
    vf(i,1)=str2num(filename(underlocs(5)+1:underlocs(6)-1));
    gammaexp(i,1)=str2num(filename((strfind(filename,'G')+1):(strfind(filename,'J')-1)));
    J(i,1)=str2num(filename((strfind(filename,'J')+1):(max(strfind(filename,'S'))-1)));
    %rep(i)=str2num(filename(max(strfind(filename,'p'))+1:underlocs(end)-1));
    rep(i,1)=str2num(filename((min(strfind(filename,'R')+1)):(strfind(filename,'Dx')-1)));  
   setnum(i,1)=str2num(filename(strfind(filename,'S')+1:strfind(filename,'t')-1));
    stepsize(i,1)=str2num(filename((strfind(filename,'Dx')+2):(strfind(filename,'G')-1)));
    t(i,1)=str2num(filename((strfind(filename,'t')+1):(max(strfind(filename,'R'))-1)));
     
end

maxNdrops=max(N);
raw=zeros(length(dirobj),maxNdrops);

for i = 1:length(fnames)
    filename=fnames{i};
    tmp=load([directory filename]);
raw(i,1:length(tmp))=tmp;
end
save([directory 'simulation_alldata_workspace.mat'])
