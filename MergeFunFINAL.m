%this function runs the simulation given a set of inputs
function MergeFunFINAL(replicates,Nstart,J,injectionstart,injectionend,T,vf,stepsize,alpha,Vmin,Vmaxstart,Dcoefexp,outdir,rep,outTs,LatticeFlag,setnum,prerunflag)
               mergeout={};
        N=(injectionend-injectionstart+1)*J+Nstart;
%draw initial Rs from a continuous uniform distribution with a minimum
        %Vmin=1;
        %Vmax=2*(Rstart.^3)-Vmin;
        Vstarts=Vmin+(Vmaxstart-Vmin).*rand(1,N); %uniform IC
%        Vstarts=10.^(log10(Vmin)+log10((Vmaxstart-Vmin)).*rand(1,N)); %log uniform IC
        Rstarts=Vstarts.^(1/3);
        L=((4/3)*pi*sum(Rstarts.^3)/vf).^(1/3); %set box size based on volume fraction and droplet number
        Rstart=mean(Rstarts);
        volbinedges=-.2:.2:6;
        
        if J>1 %J must be an integer or the reciprocal of an integer, for this implementation!
            injectiontimes=repmat((injectionstart):injectionend,1,J)';
        else 
            injectiontimes=(injectionstart:1/J:injectionend)';
        end
        injectionIDs=((Nstart+1):N)';   
        R=NaN(1,N);%preallocate radius outputs
        histbinedges=((0:1:(N+2)))*(Vmaxstart+Vmin)/2+0.5;
        Vmeanhist=zeros(1,length(histbinedges)-1);
        stepcount=0;
        R(1,1:Nstart)=Rstarts(1,1:Nstart);
        if LatticeFlag==1
            latticepointnum=ceil(N.^(1/3)); %generates lattice initial condition--if N isn't a perfect cube vf will be slightly too high. 
            [xgrid ygrid zgrid]=ndgrid(1:latticepointnum,1:latticepointnum,1:latticepointnum); 
            x0=reshape(L*(xgrid-.5)/latticepointnum, [1,N]);
            y0=reshape(L*(ygrid-.5)/latticepointnum, [1,N]);
            z0=reshape(L*(zgrid-.5)/latticepointnum, [1,N]);       
        else 
            x0= L*rand(1,N);
            y0= L*rand(1,N);
            z0= L*rand(1,N);
        end
            distvec=distcalculate(x0',y0',z0',L);
            radiimat=R(1,:)+R(1,:)';
            [p q]=find(squareform(distvec)<radiimat);
            mergepairs=unique([p(p~=q) q(p~=q)],'rows');
            while isempty(mergepairs)==0
                [p q]=find(squareform(distvec)<radiimat);
                R1=R(1,mergepairs(1,1));
                R2=R(1,mergepairs(1,2));
                R(1,mergepairs(1,2))=R1;%(nansum([R1.^3,R2.^3])).^(1/3);
                R(1,mergepairs(1,1))=NaN;
                radiimat=R(1,:)+R(1,:)';
                x0(1,mergepairs(1,2))= (R1.*x0(1,mergepairs(1,1)) +R2.*x0(1,mergepairs(1,2)))./(R1+R2);
                y0(1,mergepairs(1,2))= (R1.*y0(1,mergepairs(1,1)) +R2.*y0(1,mergepairs(1,2)))./(R1+R2);
                z0(1,mergepairs(1,2))= (R1.*z0(1,mergepairs(1,1)) +R2.*z0(1,mergepairs(1,2)))./(R1+R2);
                [p q]=find(squareform(distvec)<radiimat);
                mergepairs=unique([p(p~=q) q(p~=q)],'rows');        
            end
        %start with initial conditions     
        xtraj(1,:)=x0;
        ytraj(1,:)=y0;
        ztraj(1,:)=z0;
        %"burnin" with fast dx
        if prerunflag==1
        dxburn=.03;
Tburn=2e4+1;
alphapre=1;
            Xburn=zeros(Tburn,N);
            Yburn=zeros(Tburn,N);
            Zburn=zeros(Tburn,N);
        for i = 1:N
            Xburn(:,i)=wfbm(alphapre/2,Tburn);
            Yburn(:,i)=wfbm(alphapre/2,Tburn);
            Zburn(:,i)=wfbm(alphapre/2,Tburn);
        end   
        xscale=sqrt(3)*std(diff(Xburn,1),1);
        yscale=sqrt(3)*std(diff(Yburn,1),1);
        zscale=sqrt(3)*std(diff(Zburn,1),1);
        xburn=dxburn*Xburn./xscale;
        yburn=dxburn*Yburn./yscale;
        zburn=dxburn*Zburn./zscale;

        xburnsteps=diff(xburn,1);
        yburnsteps=diff(yburn,1);
        zburnsteps=diff(zburn,1);      
for t=1:(Tburn-1)
            %perform the injection
            R(1,injectionIDs(find(injectiontimes==t)))=Rstarts(1,injectionIDs(find(injectiontimes==t)));
            xtraj(1,injectionIDs(find(injectiontimes==t)))= L*rand(1,length(find(injectiontimes==t)));
            ytraj(1,injectionIDs(find(injectiontimes==t)))= L*rand(1,length(find(injectiontimes==t)));
            ztraj(1,injectionIDs(find(injectiontimes==t)))= L*rand(1,length(find(injectiontimes==t)));
            xtraj(1,:)=mod(xtraj(1,:)+((1./R(1,:)).^(Dcoefexp/2)).*xburnsteps(t,:),L);
            ytraj(1,:)=mod(ytraj(1,:)+((1./R(1,:)).^(Dcoefexp/2)).*yburnsteps(t,:),L);
            ztraj(1,:)=mod(ztraj(1,:)+((1./R(1,:)).^(Dcoefexp/2)).*zburnsteps(t,:),L);
            %R(t+1,:)=R(t,:);
            radiimat=R(1,:)+R(1,:)';
            distvec=distcalculate(xtraj(1,:)',ytraj(1,:)',ztraj(1,:)',L);
            [p q]=find(squareform(distvec)<radiimat);
            mergepairs=unique([p(p~=q) q(p~=q)],'rows');
            %perform the merger
            while isempty(mergepairs)==0
                [p q]=find(squareform(distvec)<radiimat);
                R1=R(1,mergepairs(1,1));
                R2=R(1,mergepairs(1,2));
                R(1,mergepairs(1,2))=R1;%(nansum([R1.^3,R2.^3])).^(1/3);%(nansum([R1.^3,R2.^3])).^(1/3);
                R(1,mergepairs(1,1))=NaN; 
                radiimat=R(1,:)+R(1,:)';
                xtraj(1,mergepairs(1,2))= (R1.*xtraj(1,mergepairs(1,1)) +R2.*xtraj(1,mergepairs(1,2)))./(R1+R2);
                ytraj(1,mergepairs(1,2))= (R1.*ytraj(1,mergepairs(1,1)) +R2.*ytraj(1,mergepairs(1,2)))./(R1+R2);
                ztraj(1,mergepairs(1,2))= (R1.*ztraj(1,mergepairs(1,1)) +R2.*ztraj(1,mergepairs(1,2)))./(R1+R2);
                distvec=distcalculate(xtraj(1,:)',ytraj(1,:)',ztraj(1,:)',L);
                [p q]=find(squareform(distvec)<radiimat);
                mergepairs=unique([p(p~=q) q(p~=q)],'rows');        
            end 
end              
        end
        %%do the actual simulation    
            X=zeros(T,N);
            Y=zeros(T,N);
            Z=zeros(T,N);
        for i = 1:N
            X(:,i)=wfbm(alpha/2,T);
            Y(:,i)=wfbm(alpha/2,T);
            Z(:,i)=wfbm(alpha/2,T);
        end   
        xscale=sqrt(3)*std(diff(X,1),1);
        yscale=sqrt(3)*std(diff(Y,1),1);
        zscale=sqrt(3)*std(diff(Z,1),1);
        x=stepsize*X./xscale;
        y=stepsize*Y./yscale;
        z=stepsize*Z./zscale;
        xsteps=diff(x,1);
        ysteps=diff(y,1);
        zsteps=diff(z,1);        
        for t=1:(T-1)
            %perform the injection
            R(1,injectionIDs(find(injectiontimes==t)))=Rstarts(1,injectionIDs(find(injectiontimes==t)));
            xtraj(1,injectionIDs(find(injectiontimes==t)))= L*rand(1,length(find(injectiontimes==t)));
            ytraj(1,injectionIDs(find(injectiontimes==t)))= L*rand(1,length(find(injectiontimes==t)));
            ztraj(1,injectionIDs(find(injectiontimes==t)))= L*rand(1,length(find(injectiontimes==t)));
            xtraj(1,:)=mod(xtraj(1,:)+((1./R(1,:)).^(Dcoefexp/2)).*xsteps(t,:),L);
            ytraj(1,:)=mod(ytraj(1,:)+((1./R(1,:)).^(Dcoefexp/2)).*ysteps(t,:),L);
            ztraj(1,:)=mod(ztraj(1,:)+((1./R(1,:)).^(Dcoefexp/2)).*zsteps(t,:),L);
            %R(t+1,:)=R(t,:);
            radiimat=R(1,:)+R(1,:)';
            distvec=distcalculate(xtraj(1,:)',ytraj(1,:)',ztraj(1,:)',L);
            [p q]=find(squareform(distvec)<radiimat);
            mergepairs=unique([p(p~=q) q(p~=q)],'rows');
            %calculate and report statistics of merging pairs...this MUST be
            %done before entering the while loop. 
            if isempty(mergepairs)==0
            vols=(R.^3);
            logvols=log10(vols);
            V1s=vols(1,mergepairs(:,1))';
            V2s=vols(1,mergepairs(:,2))';
            logV1s=log10(V1s);
            logV2s=log10(V2s);
            volbinedges=-.2:.2:3.6;
            volbincenters=(volbinedges(1:end-1)+volbinedges(2:end))/2;           
            [nstmp x]=histcounts(logvols,'BinEdges',volbinedges);
            Nclusters=sum(nstmp);
            V1bin=zeros(size(V1s,1),1);
            V2bin=zeros(size(V2s,1),1);
            for v1ind=1:length(V1s)
            V1bin(v1ind,1)=max(find(logV1s(v1ind)>volbinedges));
            end
            for v2ind=1:length(V2s)
            V2bin(v2ind,1)=max(find(logV2s(v2ind)>volbinedges));
            end
            NV1=nstmp(V1bin)';
            NV2=nstmp(V2bin)';              
            mergeout{t,1}=[(t-1)*ones(size(mergepairs,1),1) mergepairs V1s V2s NV1 NV2 Nclusters*ones(size(mergepairs,1),1) V1bin V2bin];        
            end
            %perform the merger
            while isempty(mergepairs)==0
                [p q]=find(squareform(distvec)<radiimat);
                R1=R(1,mergepairs(1,1));
                R2=R(1,mergepairs(1,2));
                R(1,mergepairs(1,2))=R1;%(nansum([R1.^3,R2.^3])).^(1/3);%(nansum([R1.^3,R2.^3])).^(1/3);
                R(1,mergepairs(1,1))=NaN; 
                radiimat=R(1,:)+R(1,:)';
                xtraj(1,mergepairs(1,2))= (R1.*xtraj(1,mergepairs(1,1)) +R2.*xtraj(1,mergepairs(1,2)))./(R1+R2);
                ytraj(1,mergepairs(1,2))= (R1.*ytraj(1,mergepairs(1,1)) +R2.*ytraj(1,mergepairs(1,2)))./(R1+R2);
                ztraj(1,mergepairs(1,2))= (R1.*ztraj(1,mergepairs(1,1)) +R2.*ztraj(1,mergepairs(1,2)))./(R1+R2);
                distvec=distcalculate(xtraj(1,:)',ytraj(1,:)',ztraj(1,:)',L);
                [p q]=find(squareform(distvec)<radiimat);
                mergepairs=unique([p(p~=q) q(p~=q)],'rows');        
            end   
            [Rhisttmp x]=histcounts(R.^3,'BinEdges',histbinedges);
                       Vmeanhist=Rhisttmp+Vmeanhist;
                       stepcount=stepcount+1;
                   if ismember(t,outTs)==1 

                         if LatticeFlag==1
                                         %suppressed outputs give
                                         %coordinates in xyz space
          %  outfilenameX=[outdir 'LATTICEalpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'X.csv'];
          %  outfilenameY=[outdir 'LATTICEalpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'Y.csv'];
          %  outfilenameZ=[outdir 'LATTICEalpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'Z.csv'];
            outfilenameR=[outdir 'LATTICEalpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'J' num2str(J) 'gamma' num2str(Dcoefexp) 'R.csv'];
     %       outfilenameRstats=[outdir 'LATTICEalpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'Rstats.csv'];
                         else       
        %outfilenameX=[outdir 'alpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' um2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'X.csv'];
    %outfilenameY=[outdir 'alpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'Y.csv'];
    %outfilenameZ=[outdir 'alpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 't' num2str(t) 'gamma' num2str(Dcoefexp) 'Z.csv'];
    outfilenameR=[outdir 'A_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_R' num2str(rep) 'Dx' num2str(stepsize) 'G' num2str(Dcoefexp) 'J' num2str(J) 'S' num2str(setnum)  't' num2str(t)   'R.csv'];
    outfilenameVmean=[outdir 'A_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_R' num2str(rep) 'Dx' num2str(stepsize) 'G' num2str(Dcoefexp) 'J' num2str(J) 'S' num2str(setnum)  't' num2str(t)   'Vmean.csv'];
    
    %outfilenameRstats=[outdir 'alpha_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_rep' num2str(rep) 'Rstats.csv'];
                         end
                xtrajout=xtraj;
                ytrajout=ytraj;
                ztrajout=ztraj;
                Rout=R;
                outputhist=Vmeanhist/stepcount;
           %   parsave(outfilenameX,xtrajout)
           %     parsave(outfilenameY,ytrajout)
           %     parsave(outfilenameZ,ztrajout)
                parsave(outfilenameR,Rout)     
               % parsave(outfilenameVmean,outputhist)
                Vmeanhist=zeros(1,length(histbinedges)-1);
           % stepcount=0;
                   end 

        end
               mergeoutmat=cell2mat(mergeout);
            outfilenameMerge=[outdir 'A_' num2str(alpha) '_N_' num2str(N) '_V_' num2str(vf) '_R' num2str(rep) 'Dx' num2str(stepsize) 'G' num2str(Dcoefexp) 'J' num2str(J) 'S' num2str(setnum) '_mergers.csv'];
                   parsave(outfilenameMerge,mergeoutmat);

    function parsave(fname, x)
      save(fname, 'x','-ascii')
    end

    function distpd =distcalculate(x0,y0,z0,L)
    [xdist]=pdist(x0);
    xdistcomplement=L-xdist;
    xdistpd=xdist;
    xdistpd(xdist>(L/2))=xdistcomplement(xdist>(L/2));


    [ydist]=pdist(y0);
    ydistcomplement=L-ydist;
    ydistpd=ydist;
    ydistpd(ydist>(L/2))=ydistcomplement(ydist>(L/2));

    [zdist]=pdist(z0);
    zdistcomplement=L-zdist;
    zdistpd=zdist;
    zdistpd(zdist>(L/2))=zdistcomplement(zdist>(L/2));

    distpd=sqrt(zdistpd.^2+ydistpd.^2+xdistpd.^2);
    end
    end