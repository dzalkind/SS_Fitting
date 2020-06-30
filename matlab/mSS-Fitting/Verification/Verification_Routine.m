%% Verification
%Routine to verify the fiting of the state-space models implementation,
%comparing the different fitting methods.
%Models available:
% 1-Spar with R^2=0.97 (use spar_097 as platform name)
% 2-Spar with R^2=0.99 (use spar_099 as platform name)
% 3-Semi with R^2=0.97 (use marin_semi_097 as platform name)
%
% Tiago Duarte
% Instituto Superio Tecnico - IST Lisbon, Portugal
%
% in collaboration with:
% WavEC Offshore Renewables (www.wavec.org)
% NREL (www.wind.nrel.gov)

clear all
name=input('Which platform results do you wish to run?','s');
addpath(name);
load(strcat(name,'_K.mat')); %Includes the vectors Ks, Kt, wi, t
clear Ks

%% Get the mass and hidrostatic matrices

% Reads mass file
fid = fopen(strcat(name,'_mass.txt'));
Rad = textscan(fid,'%s','delimiter','%','commentstyle','%','Headerlines',1);
fclose(fid);
% str2num(Rad{1}{3})str2num(Rad{1}{4})str2num(Rad{1}{5})
V=str2num(Rad{1}{6});                       %Displaced volume
CG=str2num(Rad{1}{1});                       %Center of Gravity (Xg,Yg,Zg)
RG=[1 0. 0.; ...
    0. 1 0.; ...
    0. 0. 1];  %Radii of gyration
Mass=str2num(Rad{1}{2}); 

% Get the Hydrostatic matrix
   ULEN=1;
   fid0=fopen(strcat(name,'.HST'));
   Hst=textscan(fid0,'%n %n %f','commentstyle','W');
   fclose(fid0);
   m=max(Hst{2}); frm='%3.0f\t%2.6E\t ';
   HS(m,2*m)=0; M(m,2*m)=0;   
   for s=2:m, frm=horzcat(frm,'%3.0f\t%2.6E\t'); end
   frm=horzcat(frm,'\n');
%    fid0a=fopen('HydroStatic.log','a+');
%    fid0b=fopen('MassMatrix.log','a+');
   for ij=1:1:m, HS(:,2*ij-1)=(10+ij:10:10*m+ij); M=HS; end
   for k=1:1:length(Hst{2}) 
       b=horzcat(Hst{1}(k),Hst{2}(k));
       x=rem(b(1),6); if x==0; x=6; end
       y=rem(b(2),6); if y==0; y=6; end
       switch x & y
             case 3 & 3,          DIM=ULEN^2;
             case(3 & 4)|(3 & 5), DIM=ULEN^3;
             otherwise,           DIM=ULEN^4;
       end
       HS(b(2),2*b(1))=Hst{3}(k)*1025*9.81*DIM;
       HS(b(1),2*b(2))=Hst{3}(k)*1025*9.81*DIM;
   end
   a=1025*V(1)*eye(3,3); b=RG.*abs(RG);
   c=diag([str2num(Rad{1}{3}) str2num(Rad{1}{4}) str2num(Rad{1}{5})]);%1025*V(1)*[0.0 CG(3) -CG(2); -CG(3) 0.0 CG(1); CG(2) -CG(1) 0.0];
   MM=[a,b;b,c];
   for r=1:1:6
       for s=1:1:6
           M(s,r*2)=MM(s,r); 
       end
   end

HySt = [1 zeros(1,12); HS zeros(6,1)]; 
Ine =[1 zeros(1,12); M zeros(6,1)];   

% Get the AddedMass, INFAddedMass and Damping Matrices
A=[];
B=[];
   fid1 = fopen(strcat(name,'.1'));
   OPT1 = textscan(fid1,'%f %n %n %f %f','commentstyle','W');
   fclose(fid1); 
   ctr=OPT1{1}(1); %-1 if 0 freq, 0 if Inf freq
   if ctr~=0, 
       n=1; 
       while OPT1{1}(n)==OPT1{1}(n+1), 
           n=n+1; 
       end
   end
   if ctr==0, %Inf freq
       n=length(OPT1{1}); 
   end
   m=max(OPT1{2}); frm='%3.0f\t%2.6E\t ';
   AM(m,2*m)=0; DM(m,2*m)=0;   
   for s=2:m, frm=horzcat(frm,'%3.0f\t%2.6E\t'); end
   frm=horzcat(frm,'\n');
   for ij=1:1:m, AM(:,2*ij-1)=(10+ij:10:10*m+ij); end
   if ctr<0,
       OPT1{5}(1:n)=0.0; 
   end
   if ctr==0, fid1a=fopen('INFAddedMass.log','a+'); else
%       fid1a=fopen('AddedMass.log','a+');
%       fid1b=fopen('Damping.log','a+');
      for ij=1:1:m, DP(:,2*ij-1)=(10+ij:10:10*m+ij); end
   end
   for r=0:n:length(OPT1{1})-n
       for k=1+r:1:n+r  
           b=horzcat(OPT1{2}(k),OPT1{3}(k));
           AM(b(2),2*b(1))=OPT1{4}(k);
           AM(b(1),2*b(2))=OPT1{4}(k);
           if ctr~=0.0
              DP(b(2),2*b(1))=OPT1{5}(k);
              DP(b(1),2*b(2))=OPT1{5}(k);
           end
           d=mod(k,m);
           if d==0
              switch ctr & r
                  %ver melhor a freq inf & zero
                  case (-1 & 0), FREQ=0.0; frms=' %3.3f\t%3.3f\n';
                  case (0 & 0), FREQ='INF FREQ'; frms=' %8s\t%3.3f\n';
                  otherwise, FREQ=2*pi/OPT1{1}(k); frms=' %3.3f\t%3.3f\n';
              end
           end
       end

       if isfinite(FREQ)
            A=[A; FREQ ULEN zeros(1,11); AM zeros(6,1)];
            B=[B; FREQ ULEN zeros(1,11); DP zeros(6,1)];
       else
            Ainf=AM;
       end
   end    

gp.dof=[1,2,3,4,5,6];

%%%%%%%%%%%%%%%%%%%%%% MASS and HYDROSTATIC ARRAYS %%%%%%%%%%%%%%%%%%%%%%%%
RMass=Inertia(gp.dof,Ine,Ainf);     %Real mass = mass + inf. added mass
Hst=HydroStatic(gp.dof,HySt);       %Hydrostatic coefficients array
% 
HM=[diag(HySt(2:end,2:2:end)) diag(Ine(2:end,2:2:end))];

clear AM DP
%Radiation Matrices
dof=[1,2,3,4,5,6];
ldof=6;
ii=ldof;
kmax=length(A(:,1))/(ldof+1);%number of frequency
AM(kmax,ii,ii)=0; %Initialize matrix
AMinf=AM; %Initialize infinite Added Mass
D=AM;%Initialize matrix
K=AM;%Initialize matrix
Dd=AM;%Initialize matrix
AMd=AM;%Initialize matrix
AMinfd=zeros(6,6);%Initialize matrix
w=zeros(1,kmax);%Initialize matrix

for l=1:1:ii
    kl=dof(l);
    for j=1:1:ii
        kj=dof(j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: For more than 6DoF is advisable  %
        % to check the dimensionless method      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x=rem(l,6); if x==0; x=6; end
        y=rem(j,6); if y==0; y=6; end
        switch (x & y)
               case (x<=3 & y<=3), n=3;
               case (x<=3 & y>3)|(x>3 & y<=3), n=4;
               case (x>3 & y>3), n=5;
        end
        for k=1:kmax % for each frequency
            w(k)=A((k-1)*(ldof+1)+1,1);
            DimAM=1000*(A(1,2))^n;
            DimD=1000*(A(1,2))^n*w(k);
            AMinf(l,j)=(Ainf(kl,2*kj));
            AM(k,l,j)=(A((kl+1)+((k-1)*(ldof+1)),2*kj));
            D(k,l,j)=(B((kl+1)+((k-1)*(ldof+1)),2*kj));
            K(k,l,j)=(D(k,l,j)*DimD)+1i*w(k)*(AM(k,l,j)-AMinf(l,j))*DimAM;
            Dd(k,l,j)=D(k,l,j)*DimD;
            AMd(k,l,j)=AM(k,l,j)*DimAM;
            AMinfd(l,j)=AMinf(l,j)*DimAM;
        end
    end
end
%clear A Ainf B AM AMinf D   
%clear x y l kl j kj n DimAM DimD nDoF 
expnt=8; %Interpolated sample number {2^expnt}                         
wi = 0.0:w(k)/(2^expnt-1):w(k); 
np=length(wi);
Ks(np,ii,ii)=0;
Ks(:,:,:)= interp1(w,K,wi,'linear','extrap');
%DKs=zeros(1);
Ds(:,:,:)= interp1(w,Dd,wi,'linear','extrap');
AMds(:,:,:)= interp1(w,AMd,wi,'linear','extrap');

DKs= reshape(max(D),ii,ii);
DKs=max(diag(DKs));

clear A B C x
for s=1:1:ii
    for r=s:1:ii 
        if r~=s
           if max(abs(D(:,s,r)))<=0.0005*DKs,
              %if max(abs(Ks(:,s,r)))<=0.05*DKs, 
              Ks(:,s,r)=0.0; 
           end
        else
           if max(abs(D(:,s,s)))<=0.0001*DKs,
              %if max(abs(Ks(:,s,s)))<=0.0001*DKs, 
              Ks(:,s,s)=0.0; 
%               disp([identify(dof(s)),' - Unimportant DoF']);
           end 
        end
    end
end
%% Load the SS Models
model{1}=[name '_FREQ.ss'];
model{3}=[name '_Time.ss'];
model{2}=[name '_FDI.ss'];
model{4}=[name '_Time2.ss'];

for ii=1:4
    
    % Get up the model matrices
    Matrices{ii}=importdata(model{ii});
    
    ns=size(Matrices{ii},2);%Number of radiation states
    
    A{ii}=Matrices{ii}(1:ns,:); %A Matrix
    B{ii}=Matrices{ii}(ns+1:ns+ns,1:6); %B Matrix
    
    if sum(sum(isnan(B{ii}(:))))>0%Missing one dof, assuming it is Yaw
        B{ii}(isnan(B{ii}))=0;
        C{ii}=Matrices{ii}(end-4:end,:);
        C{ii}(end+1,:)=0;%C Matrix
    else
    C{ii}=Matrices{ii}(end-5:end,:); %C Matrix
    end
    D=[];
    
end

%% Reference K(w),K(t), Hf2d and Hf2v
COLOR={'b','g','r','c'};

%Get the non zero entries of K
q=0;
ij=[];
for s=1:1:6 % for each...
    for r=s:1:6 % entry of the matrix ..
        if max(abs(Ks(:,s,r)))>0.0
            ij=[ij; s,r];
            q=q+1;            
        end        
    end
end            

%Reference plots
% Ks
scrsz = get(0,'ScreenSize');
h=figure('name','Transfer Function','OuterPosition',[1 1 scrsz(3) scrsz(4)]);
for ii=1:q
    
    hh=subplot(floor(q/2)+1,2,ii);
	plot(wi,abs(Ks(:,ij(ii,1),ij(ii,2))),'k','LineWidth',2)
    hold on
%     if ii<q
%         set(hh,'XTick',[]);
%     else
        xlabel('Freq [rad/s]');
%     end
    
    ylabel(strcat('K_{',num2str(ij(ii,1)),num2str(ij(ii,2)),'}(jw)'))
end
            
%Kt
k=figure('name','Impulse response Function','OuterPosition',[1 1 scrsz(3) scrsz(4)]);
for ii=1:q
    
    kk=subplot(floor(q/2)+1,2,ii);
	plot(stime,Kt(:,ii),'k','LineWidth',2)
    hold on
%     if ii<q
%         set(kk,'XTick',[]);
%     else
        xlabel('Time [s]')
%     end
    ylabel(strcat('K_{',num2str(ij(ii,1)),num2str(ij(ii,2)),'}(t)'))
end

%Hf2d & Hf2v
for ww=2:length(wi)%for each frequency
  
    AAA(:,:)=AMds(ww,:,:);
    DDD(:,:)=Ds(ww,:,:);
    
    G=((1i*wi(ww))^2.*eye(6,6)+(RMass+AMinfd)^-1*Hst)^-1*(RMass+AMinfd)^-1*1i*wi(ww);

    H(:,:)=Ks(ww,:,:);        
%     Hf2d(:,:,ww)=(1i*wi(ww)).^-1.*(eye(6,6)+G*H(:,:))^-1*H(:,:);
%     Hf2v(:,:,ww)=(eye(6,6)+G*H(:,:))^-1*H(:,:);

    
    
    Hf2d(:,:,ww)=(-wi(ww)^2.*(RMass+AAA)-1i*wi(ww).*DDD+Hst)^-1;
    Hf2v(:,:,ww)=1i*wi(ww)*Hf2d(:,:,ww);
end


cc=figure('name','Added Mass','OuterPosition',[1 1 scrsz(3) scrsz(4)]);
for ii=1:q
    
    hh=subplot(floor(q/2)+1,2,ii);
    g(:,1)=AMds(:,ij(ii,1),ij(ii,2));
	plot(wi(2:end),g(2:end),'k','LineWidth',2)
    hold on
%     if ii<q
%         set(hh,'XTick',[]);
%     else
        xlabel('Freq [rad/s]');
%     end
    
    ylabel(strcat('AM_{',num2str(ij(ii,1)),num2str(ij(ii,2)),'}(jw)'))
end

dd=figure('name','Damping','OuterPosition',[1 1 scrsz(3) scrsz(4)]);
for ii=1:q
    
    hh=subplot(floor(q/2)+1,2,ii);
    g(:,1)=Ds(:,ij(ii,1),ij(ii,2));
	plot(wi,g,'k','LineWidth',2)
    hold on
%     if ii<q
%         set(hh,'XTick',[]);
%     else
        xlabel('Freq [rad/s]');
%     end
    
    ylabel(strcat('B_{',num2str(ij(ii,1)),num2str(ij(ii,2)),'}(jw)'))
end

%% Comparison of the results

for ii=1:4 %For each model
    
    %Create the global SS system
    sys=ss(A{ii},B{ii},-C{ii},0,0);
    
    %General calculations
    [mag,phase] = bode(sys,wi); %Frequency response (tf)
    [y,t] = impulse(sys,stime); %Impulse response function
    [H] = freqresp(sys,wi);
      
    
    for ww=2:length(wi)%for each frequency

        G=((1i*wi(ww))^2.*eye(6,6)+(RMass+AMinfd)^-1*Hst)^-1*(RMass+AMinfd)^-1*1i*wi(ww);
        
%         Hf2de(:,:,ww)=(1i*wi(ww)).^-1.*(eye(6,6)+G*H(:,:,ww))^-1*H(:,:,ww);
%         Hf2ve(:,:,ww)=(eye(6,6)+G*H(:,:,ww))^-1*H(:,:,ww);
        ADDEDM(:,:,ww)=imag(wi(ww)^-1*H(:,:,ww))+AMinfd;
        DAMP(:,:,ww)=real(H(:,:,ww));
        
        Hf2de(:,:,ww)=(-wi(ww)^2.*(RMass+ADDEDM(:,:,ww))-1i*wi(ww).*DAMP(:,:,ww)+Hst)^-1;
        Hf2ve(:,:,ww)=1i*wi(ww)*Hf2d(:,:,ww);

    end
                      
    %Ks Comparisons
    figure(h)
    for jj=1:q
        hh=subplot(floor(q/2)+1,2,jj);
        g(:,1)=mag(ij(jj,1),ij(jj,2),:);
        plot(wi',g,COLOR{ii})
        hold on
        if jj==q && ii==4
                legend('Orig','Freq. Domain', 'FDI','TD - Prony','TD - IMP2SS','Location','EastOutside','Orientation','horizontal')
        end
        
        %Get the values of Rsquare for each Kij
        res=abs(Ks(:,ij(jj,1),ij(jj,2)))-g(:,1);                     %Residuals
        avrg=mean(abs(Ks(:,ij(jj,1),ij(jj,2))));             %Sum of squares about the mean
        sst=sum(((abs(Ks(:,ij(jj,1),ij(jj,2)))-avrg).*conj(abs(Ks(:,ij(jj,1),ij(jj,2)))-avrg)));  %Weighted average value
        Rsquare_freq(ii,jj)=1-sum(res.^2)/sst; %Sum of squares about the mean
             
    end  
    
    %Added mass comparisons
    figure(cc)
    for jj=1:q
        hh=subplot(floor(q/2)+1,2,jj);
        g(:,1)=ADDEDM(ij(jj,1),ij(jj,2),:);
        plot(wi(2:end),(g(2:end)),COLOR{ii})
        hold on
        if jj==q && ii==4
                legend('Orig','Freq. Domain', 'FDI','TD - Prony','TD - IMP2SS','Location','EastOutside','Orientation','horizontal')
        end
        

    end  

        %Damping comparisons
    figure(dd)
    for jj=1:q
        hh=subplot(floor(q/2)+1,2,jj);
        g(:,1)=real(H(ij(jj,1),ij(jj,2),:));
        plot(wi',(g),COLOR{ii})
        hold on
        if jj==q && ii==4
                legend('Orig','Freq. Domain', 'FDI','TD - Prony','TD - IMP2SS','Location','EastOutside','Orientation','horizontal')
        end
        
             
    end  
    
    
    
    %Kt Comparisons
    figure(k)
    for jj=1:q
        hh=subplot(floor(q/2)+1,2,jj);
        plot(stime,y(:,ij(jj,1),ij(jj,2)),COLOR{ii})
        hold on
        axis([0 20 -inf +inf])
        if jj==q && ii==4
                legend('Orig','Freq. Domain', 'FDI','TD - Prony','TD - IMP2SS','Location','EastOutside','Orientation','horizontal')
        end
        
        %Get the values of Rsquare for each Kij
        res=abs(Kt(:,jj)-y(:,ij(jj,1),ij(jj,2)));                     %Residuals
        avrg=mean(Kt(:,jj));             %Sum of squares about the mean
        sst=sum((Kt(:,jj)-avrg).*conj(Kt(:,jj)-avrg));  %Weighted average value
        Rsquare_imp(ii,jj)=1-sum(res.^2)/sst; %Sum of squares about the mean
        
    end
  
   
    %Get the number of states per Kij
    x{ii}=[]; 
    qq=1;
    x{ii}(1,q)=0;
    for jj=1:length(C{ii})
        
        if jj==1;
            x{ii}(1,qq)=1;
            
            
        elseif qq~=q
            
            if ij(qq+1,1)==ij(qq+1,2)
                    
                if C{ii}(ij(qq+1,1),jj)==0  
                    x{ii}(1,qq)=x{ii}(1,qq)+1;
                else
                    qq=qq+1;
                    x{ii}(1,qq)=x{ii}(1,qq)+1;
                end
                
            elseif  ij(qq+1,1)~=ij(qq+1,2)
                
                if (C{ii}(ij(qq+1,2),jj)==0)
                    x{ii}(1,qq)=x{ii}(1,qq)+1;
                else
                    qq=qq+1;
                    x{ii}(1,qq)=x{ii}(1,qq)+1;
                end
            end      
        else
            x{ii}(1,qq)=x{ii}(1,qq)+1;
        end
    end
    for jj=1:q
        if ij(jj,1)~=ij(jj,2)
            x{ii}(1,jj)=x{ii}(1,jj)/2;
        end
    end
    clear y mag
    
end

MARKER={'.','.','.','.'};

figure('name','R^2 - Frequency response')
for ii=1:4
    plot(x{ii},Rsquare_freq(ii,:),strcat(MARKER{ii},COLOR{ii}),'MarkerSize',20)
    hold on    
    
end
axis([1 max([x{1},x{2},x{3},x{4}])+1 -inf 1])
grid on
legend(['Freq. Domain - ', num2str(length(A{1})),' states'],['FDI - ', num2str(length(A{2})),' states'],['TD - Prony - ', num2str(length(A{3})),' states'], ['TD - IMP2SS - ', num2str(length(A{4})),' states'],'Location','SouthEast')
xlabel('Number of states per K_{ij}')
ylabel('R^2')
title('R^2 - Frequency response')

figure('name','R^2 - Impulse response Function')
for ii=1:4
    plot(x{ii},Rsquare_imp(ii,:),strcat(MARKER{ii},COLOR{ii}),'MarkerSize',20)
    hold on
%     axis([1 6 -inf 1])
end
axis([1 max([x{1},x{2},x{3},x{4}])+1 -inf 1])
grid on
legend(['Freq. Domain - ', num2str(length(A{1})),' states'],['FDI - ', num2str(length(A{2})),' states'],['TD - Prony - ', num2str(length(A{3})),' states'], ['TD - IMP2SS - ', num2str(length(A{4})),' states'],'Location','SouthEast')
xlabel('Number of states per K_{ij}')
ylabel('R^2')
title('R^2 - Impulse response Function')

rmpath(name);
