function [condlist Nsphere S VarS Srep VarSrep nsavg nssem ccdfavg ccdfsem Vdmean Vdvar volbincenters totvol monnum monnumsem Vdmeansem Vdvarsem] = workspace_parse_function(directory)
v0=1.5;
filename=[directory '\simulation_alldata_workspace.mat'];
load(filename)
condmat = table(alphaexp,N,vf,t,gammaexp,J,stepsize);
condlist=unique(condmat,'rows');
logvolbinedges=.1:0.2:3;%log10(0.5:100:1500.5);
%logvolbinedges=log10(0.5:100:1500.5);
logvolbincenters=log10((10.^logvolbinedges(1:end-1)+10.^logvolbinedges(2:end))/2);           
for condind=1:size(condlist,1)
    [q idx]=ismember(condmat,condlist(condind,:),'rows');
    replist=find(q==1);
    numreps(condind)=length(replist);
    tmp=raw(replist,:).^3;
    Ndrops=condlist.N(condind,:);
    volbinedges=((0:1:(Ndrops+1)))*v0+0.5;
    volbincenters=(volbinedges(1:end-1)+volbinedges(2:end))/2;
Nt(condind,1)=sum(mean(~isnan(tmp),1));
for repnum=1:size(tmp,1)
    [nstmp x]=histcounts(tmp(repnum,:),'BinEdges',volbinedges);
    avgsize(repnum)=sum(nstmp.*volbincenters.^2)./sum(nstmp.*volbincenters);
    ns{condind}(repnum,:)=nstmp;
    secondmoment(repnum)=sum(nstmp.*volbincenters.^3)./sum(nstmp.*volbincenters);
    totvoltmp(repnum)=sum(tmp(repnum,:));
    minsizetmp(repnum)=min(tmp(repnum,:));
    nsf{condind}(repnum,:)=cumsum(ns{condind}(repnum,:))./sum(ns{condind}(repnum,:));
    ccdf{condind}(repnum,:)=1-nsf{condind}(repnum,:);
    Srep(repnum,condind)=avgsize(repnum);
    VarSrep(repnum,condind)=secondmoment(repnum)-avgsize(repnum).^2;
    Nsphere(repnum,condind)=sum(~isnan(tmp(1,:)));
    
end
minsize(condind)=min(minsizetmp);
    S{condind}=nanmean(avgsize);
    S2{condind}=nanmean(secondmoment);
    VarS{condind}=nanmean(secondmoment-avgsize.^2);
    nsavg{condind}=mean([ns{condind}],1);
    nssem{condind}=sem([ns{condind}],1);
    nsfavg{condind}=mean([nsf{condind}],1);
    nsfsem{condind}=sem([nsf{condind}],1);
    ccdfavg{condind}=mean([ccdf{condind}]);
    ccdfsem{condind}=sem([ccdf{condind}],1);
    monnum(condind)=mean(sum(tmp<2,2));
    monnumsem(condind)=sem(sum(tmp<2,2),1);
    
    NS=[nsavg{condind}];
    NScell{condind}=[nsavg{condind}];
  Vdmean(condind,1)=mean(nanmean(tmp,2));
    Vdmeansem(condind,1)=sqrt(var(nanmean(tmp,2))/size(tmp,1));
  nonmontmp=tmp;
  nonmontmp(nonmontmp<1.1)=NaN;
  Vdmeannonmon(condind,1)=nanmean(nanmean(nonmontmp,2));
    Vdvar(condind,1)=mean(nanvar(tmp,1,2));
    Vdvarsem(condind,1)=sqrt(var(nanvar(tmp,1,2))/size(tmp,1));
    
totvol(condind,1)=mean(nansum(tmp,2));
    NSF{condind}=cumsum(NS)./sum(NS);
    end
Rs=reshape(tmp,[1 prod(size(tmp))]);
end
