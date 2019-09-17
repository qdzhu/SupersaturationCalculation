%*************************************************************%
% Program Name:  SP_UZ0.m                                     %
% CALCULATION OF DROPLET SIZE DISTRIBUTION IN A VERTICALLY    %
% MOVING ADIABATIC PARCEL                                     %
%                                                             %
% Atmospheric Environment Service, Toronto, 1994  (version 3) %
% Environment Canada - October 2013  (version 4)              %
% by Alexei Korolev                                           %
% email: alexei.korolev@canada.ca                             %
% tel: (416) 739-5716                                         %
%*************************************************************%
clear variables 
close all

% dry CCN distribution
k=2; %CCN type flag: 1 marine; 2 clean continental; 3 background; 4 urban
nchan=300;  %number of CCN bins
rmin=0.005; %min CCN size (um)
rmax=1;     %max CCN size (um)
% 
% INPUT VALUES 
Sup=-1;        % initial supersaturation (%) 
Tco=0;         % initial temperature (%) 
Ho=500;        % initial height (m) 
Hmax=700;      % max altitude (m)
% WARNING: for Uz>5m/s reduce dtd1,2,3,4
Vel=1;         % vertical velocity (m/s) 
%----------------------------------------------------------%
% WARNING:                                                 %
% too large dt may lead to errors in DSD calculations      %
% should be dt~<dh/Vel;  dh=0.5m                           %
% dt should be multiple of dtd1, dtd2, dtd3, dtd4          %  
%----------------------------------------------------------%
dt=.5;        % time step of external loop (s) 
tmax=(Hmax-Ho)/Vel; % duration of simulation (s) 
ro=.75e-6;    % boundary droplet size for calculation of conc.in range r>ro 
%
dtd1=.05;   Nmx1=dt/dtd1; %time step and number of cycles for droplet calc.
dtd2=.01;   Nmx2=dt/dtd2; %time step and number of cycles for droplet calc.
dtd3=.0002; Nmx3=dt/dtd3; %time step and number of cycles for droplet calc.
dtd4=.0001; Nmx4=dt/dtd4; %time step and number of cycles for droplet calc.
%
% 
%----------------------------------------------------%
% definition of coefficients used in the following   %
% calculations                                       %
[a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,...
    Mms,alpha,w,Po,To,g,k_mu,k_ml]=Constant;         %
%----------------------------------------------------%
% START VALUES 
Tc=Tco;             % initial temperature
Es=(a0+Tc*(a1+Tc*(a2+Tc*(a3+Tc*(a4+Tc*(a5+Tc*a6))))))*100;  %vapor pressure over flat water surface 
E=Es*(1+Sup*1e-2);  % initial water vapour pressure 
S=E/Es;             % initial supersaturation
T=Tc+To;            % initial temperature in Kelvin
P=Po*exp(-g*Ho/(Ra*T)); % initial pressure
P0=P;               % for file recording
Ma=1;               % weight of dry air =const
Mv=Ma*Mmv*E/(Mma*(P-E)); 
Mt=Ma+Mv;           % total weight of adiabatic bolb
dM=0; 
Rt=Rg*(Mv*Mma+Ma*Mmv)/((Ma+Mv)*Mma*Mmv);%specific gas constant of moist air
Cpt=(Ma*Cpa+Mv*Cpv)/(Ma+Mv);           % thermal capacity of moistered air
Vol=(Mv/Mmv+Ma/Mma)*Rg*T/P;            % initial volume of adiabatic parcel
Nstep=round(tmax/dt);                  % the whole cycle
% 
% CONSTANTS TEMPERATURE & PRESSURE DEPENDENT 
D=(2.26e-5+1.5e-7*Tc)*Po/P;  % diffusion coeff.
L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
K=2.424e-2+7.95e-5*Tc;       % thermal conductivity
sigma=7.564e-2-1.43e-4*Tc;   % water surface tension
kappa=K*Ra*T/(Cpa*P); 
mu=k_mu*sqrt(T*Rv); 
%
%------------------------------------------------------------------------% 
% Calculation of initial spectrum of wetted CCN  at supersaturation S(0) % 
%Init_SpCCN
[r,req,Nconc]=Init_SpCCN(S,pl,ps,Mmv,Mms,Rv,T,sigma,Vol,k,rmin,rmax,nchan);
%------------------------------------------------------------------------%


Nchan=size(r,2);                % numer of the size bins
Mlo=k_ml*sum(Nconc.*(r.^3));    % initial LWC
Mtot=Mv+Mlo;       % initial weight of total water per 1 kg of dry air

Cnc=Nconc;                  % auxillary array 'cnc'
Cnc(:,Nchan)=[];            % auxillary array 'cnc'
dr=nan*zeros(1,Nchan);      % pre-allocating memory for droplt radii change
Ss=nan*zeros(1,Nstep);      % pre-allocating memory for supesaturation
H=zeros(1,Nstep);           % pre-allocating memory for height
time=zeros(1,Nstep);        % pre-allocating memory for time
Uz=zeros(1,Nstep);          % pre-allocating memory for vertical velocity
Rmax=zeros(1,Nstep);        % pre-allocating memory for max radius
Conc=nan*zeros(1,Nstep);    % pre-allocating memory for concentration

%======================================%
%             GRAPHICS                 %
% Open figure for droplet size distrib.%
%======================================%
figure(4)
h1=axes('Units','normalized','Position',[.1 ,.12,.3,.8],...
'xlim',[Sup .5],'ylim',[Ho,Hmax],'yscale','linear','xscale','linear',...
'Box','on','Tag','Axis1');
xlabel('Supersaturation (%)'); ylabel('Height (m)'); hold on

h2=axes('Units','normalized','Position',[.51 ,.3,.47,.5],...
'Box','on','Tag','Axis2');
%'xlim',[0,15],'ylim',[0,110],'yscale','log','xscale','log',...

xlabel('Radius (um)'); ylabel('F(r) (cm-3/um)'); hold on


%******************************************% 
% SOLUTION OF DIFFERENTIAL EQUATIONS       %
%        ( main cycle )                    %
%******************************************%
for ii=1:Nstep 
    disp(ii)
    
  P=P-P*g*Vel*dt/(Rt*T); 
  T=T+(-g*Vel*dt+L*dM/Mt)/Cpt; 
  Tc=T-To;  
  E=P*Mv*Mma/(Mv*Mma+Ma*Mmv); 
  Es=(a0+Tc*(a1+Tc*(a2+Tc*(a3+Tc*(a4+Tc*(a5+Tc*a6))))))*100; 
  S=E/Es; 
  % 
  D=(2.26e-5+1.5e-7*Tc)*Po/P; 
  L=2.495e6-2.3e3*Tc; 
  K=2.424e-2+7.95e-5*Tc; 
  sigma=7.564e-2-1.43e-4*Tc; 
  kappa=K*Ra*T/(Cpa*P); 
  mu=k_mu*sqrt(T*Rv); 
  % 
  Pv=E/(Rv*T); 
  k2=(L/(Rv*T)-1)*L*Pv*D/(T*K); 
  k1=D*Pv/(pl*(k2+1)); 
  ksi=(D/(alpha*mu)+k2*kappa/(w*mu))/(k2+1); 
  B=2*sigma/(Rv*T*pl); 
  % 
  for jj=1:Nchan 
    if r(jj)>1e-6,                      dtd=dtd1; Nmx=Nmx1; 
      elseif r(jj)<=1e-6 && r(jj)>3e-7, dtd=dtd2; Nmx=Nmx2; 
      elseif r(jj)<=3e-7 && r(jj)>1e-8, dtd=dtd3; Nmx=Nmx3; 
      elseif r(jj)<=1e-8,               dtd=dtd4; Nmx=Nmx4; 
    end 
   
      for iii=1:Nmx 
        dr(jj)=dtd*k1*(S-(1-1.215/(r(jj)^3*req(1,jj)-1.815))*exp(B/r(jj)))/(r(jj)+ksi);
        r(jj)=r(jj)+dr(jj); 
      end 
  end 
 
  
  %coefficient and mass balance 
  Volrev=P/((Mv/Mmv+Ma/Mma)*Rg*T)*1e-6;  % reverse volume [cm^-3] 
  Ml=k_ml*sum(Nconc.*(r.^3)); 
  dM=Ml-Mlo; 
  Mv=Mtot-Ml; 
  Mlo=Ml; 
  Mt=Mv+Ma; 
  Rt=Rg*(Mv*Mma+Ma*Mmv)/(Mt*Mma*Mmv); 
  Cpt=(Ma*Cpa+Mv*Cpv)/Mt; 
  
  %------droplet spectrum calculation--------%
  dre=(r(2:end)-r(1:end-1))*1e6;
  jr=sum(ro>r); 
  delconc=(r(jr+1)-ro)*Cnc(jr)*1e6/dre(jr); 
  Nconcgr=Cnc./dre*1e-6;    %conc in cm-3
  Concmax=1.2*max(Nconcgr.*(ro<r(1:end-1))); 

  %-----output values------%
  Conc(ii)=(sum(Nconc.*(ro<r))+delconc)*Volrev;	%droplet conc (m-3)
  Ss(ii)=(S-1)*100;                             %supersaturation  
  time(ii)=ii*dt;                               %time
  H(ii)=Ho+time(ii)*Vel;                        %altitude 
  Uz(ii)=Vel;                                   %vertical velocity
  Rmax(ii)=r(end);

  %*******************************%
  %     droplet spectrum plot     %
  %*******************************%
  figure(4);     
  stairs(h2,r(1:end-1)*1e6,Nconcgr);  hold on
  set(gca,'Xlim',[0,10],'Ylim',[0,Concmax],'Yscale','linear','Xscale','linear'); 
  %set(gca,'Xlim',[0,12],'Ylim',[0,Concmax],'Yscale','log','Xscale','log'); 
  ylabel('F(r) (cm-3/um)'), xlabel('radius (um)'),
  hold off
  title(['time=',num2str(round(time(ii))),'s, Uz=',num2str(Vel,'%2.1f'),...
      'm/s, N=',num2str(Conc(ii),'%4.1f'),'cm^{-3}']) 
  
  if ii>1
      plot(h1,[Ss(ii-1),Ss(ii)],[H(ii-1),H(ii)],'b');
      delete(h3); %delete previous red dot
  end
  h3=plot(h1,Ss(ii),H(ii),'.r','MarkerSize',12); %current red dot
  title(h1,['T=',num2str(Tc,'%4.2f'),'C, H=',num2str(H(ii),'%5.1f'),'m, S=',num2str(Ss(ii),'%2.3f'),'%']) 
  
pause(.1)

end 
%

  figure(5)
  %----supersaturation plot----%
  subplot(121)
  plot(Ss(1:1:end),H(1:1:end),'b');  hold on
  set(gca,'ylim',[H(1),max(H)]), 
  xlabel('Supersaturation (%)'), ylabel('Height (m)') 
  title(['Start values: To=',num2str(Tco),'C, Ho=',num2str(Ho),'m, S=',num2str(Sup),'%']) 
  grid
  
  %----droplet number concentration plot-----%
  subplot(122)
  plot(Conc(1:1:end),H(1:1:end),'r');  hold on
  set(gca,'ylim',[H(1),max(H)]), 
  xlabel('Conc (cm-3)') 
  grid
  
  %----max droplet size-----%
%   subplot(133)
%   plot(Rmax*1e6,H,'g');  hold on
%   set(gca,'ylim',[H(1),max(H)]), 
%   xlabel('Rmax (um)') 

disp([num2str(Tco),'    ',num2str(max(Ss)),'    ', num2str(max(Conc))])
%
%%
function [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=Constant 
% CONSTANTS 
a0=6.107799961;     a1=4.436518521e-1; a2=1.42894580e-2; a3=2.65064847e-4; 
a4=3.031240396e-6;  a5=2.034080948e-8; a6=6.136820929e-11; 
Rg=8.317;      % universal gas constant [j/mol*k] 
Cpa=1005;      % Thermocapacity of dry air under constant pressure
Cva=718;       % Thermocapacity of dry air under constant volume
Mma=.02896;    % Molecular weight of dry air
Mmv=.01806;    % Molecular weight of water vapour
Ra=Rg/Mma;     % Specific gas constant of dry air
Rv=Rg/Mmv;     % Specific gas constant of water vapour
Cpv=1850;      % Thermocapacity of water vapour under constant pressure
pl=1000;       % Weight density of water
ps=2500;       % Weight density of CCN
Mms=.079;      % Molecular weight of CCN
alpha=1;       % Coefficient of condensation
%alpha=.03;    % Coefficient of condensation
w=1;           % Coefficient of thermal accomodation
Po=1e5;        % Atmosphere pressure (N/m2)
To=273.15;     % 0C in Kelvin deg.
g=9.81;        % acceleration of gravity
k_mu=sqrt(1/(2*pi)); 
k_ml=4*pi*pl/3; 
end

%%
%*************************************************************************%
% Calculation of initial CCN distribution for four different types of CCN:% 
% maritime, clean continental, background and urban                       %
% calculation of equilibrium droplet size distribution at initial         % 
% supersaturation S(0)<0, i.e. wetted CCN                                 %
%*************************************************************************%

function [r,req,Nconc]=Init_SpCCN(S,pl,ps,Mmv,Mms,Rv,T,sigma,Vol,k,rmin,rmax,nchan) 
% k=1;  % CCN type: maritime
% k=2;  % CCN type: clean continental
% k=3;  % CCN type: background
% k=4;  % CCN type: urban


%======================================================================
%  Mark Pinsky original code
%======================================================================
dr     = 0.001;
rn     = 0.001:dr:50;
n      = length(rn);

% coefficient for the dry CCN distribution
AA   = [340, 1000, 6400, 106000;...
         60,  800, 2300,  32000;...;
        3.1, 0.72,  3.2, 5.4];
    
BB   = [0.005, 0.008, 0.008, 0.007;...
        0.035, 0.034, 0.038, 0.027;...;
        0.31,   0.46,  0.51,  0.43];
    
CC   = [1.6, 1.6,  1.7, 1.8;...
        2.0, 2.1,  2.0, 2.16;...;
        2.7, 2.2, 2.16, 2.21];       

%for k=3 
   sp = (AA(1,k)/(sqrt(2*pi)*log10(CC(1,k))*log(10))*...
         exp(-(log10(rn/BB(1,k))).^2/(2*log10(CC(1,k))^2)))./rn+...
        (AA(2,k)/(sqrt(2*pi)*log10(CC(2,k))*log(10))*...
         exp(-(log10(rn/BB(2,k))).^2/(2*log10(CC(2,k))^2)))./rn+...
        (AA(3,k)/(sqrt(2*pi)*log10(CC(3,k))*log(10))*...
         exp(-(log10(rn/BB(3,k))).^2/(2*log10(CC(3,k))^2)))./rn;
     
%end
disp(sum(sp)*dr)    
%======================================================================
%======================================================================


%==========================================%
% Calculation of dry CCN distribution      %
%==========================================%
rccn=exp(log(rmin) :(log(rmax)-log(rmin))/nchan : log(rmax));

   sp1 = (AA(1,k)/(sqrt(2*pi)*log10(CC(1,k))*log(10))*...
         exp(-(log10(rccn/BB(1,k))).^2/(2*log10(CC(1,k))^2)))./rccn+...
        (AA(2,k)/(sqrt(2*pi)*log10(CC(2,k))*log(10))*...
         exp(-(log10(rccn/BB(2,k))).^2/(2*log10(CC(2,k))^2)))./rccn+...
        (AA(3,k)/(sqrt(2*pi)*log10(CC(3,k))*log(10))*...
         exp(-(log10(rccn/BB(3,k))).^2/(2*log10(CC(3,k))^2)))./rccn;
drccn=rccn(2:end)-rccn(1:end-1);
nccn=(sp1(1:end-1)+sp1(2:end))/2;     %middle of the size bins
nchan=size(rccn,2); 
Nconc=nccn.*drccn*Vol*1e6; % distribution of CCNs numbers in the volume Vol 

%---------------%
%  GRAPHICS     %
%---------------%
if k==1
    titl='marine CCN';
elseif k==2
    titl='clean continental CCN';
elseif k==3
    titl='background CCN';
elseif k==4
    titl='urban CCN';
end

figure(1)
loglog(rn,sp), hold on
%loglog(rccn,sp1,clr(2)), hold on
loglog(rccn,sp1,'.r')
xlabel('dry CCN radius/droplet radius (um)') 
ylabel('CCN radii distribution (cm^{-3} um^{-1})') 
title(titl)
legend('dry CCN radii distribution','truncated CCN radii distr.')
CCN_CONC=sum(nccn.*drccn); 
disp(['CCN concentration =', num2str(CCN_CONC),'cm-3'])


%-----------------------------% 
%   Wetted equilibrium CCN    %
%-----------------------------%
rccn=rccn*1e-6; %rccn convertion in meters
err=rccn*1e-3; 
del_r=rccn*.5; 
b=2*sigma/(Rv*T*pl); 
req=(pl*Mms/(2*ps*Mmv))*rccn.^(-3); 
r=((1.815)^(1/3))*req.^(-1/3)+err; 
y1=S-(1-1.215./(r.^3.*req-1.815)).*exp(b./r); 
% 
% Loop for equation roots (i.e. equilibrium radii)
for i=1:nchan   
  
  r0=r(i); y10=y1(i); y20=y1(i); 
  while y10*y20>0 
    r0=r0+del_r(i); 
    y20=y10; 
    y10=S-(1-1.215/(r0^3*req(i)-1.815))*exp(b/r0); 
  end 

  r1=r0; r2=r0-del_r(i); 
  while abs(r1-r2)>err(i) 
    r3=(r1+r2)*.5; 
    y30=S-(1-1.215/(r3^3*req(i)-1.815))*exp(b/r3); 
    if y10*y30<0,  y20=y30; r2=r3; end 
    if y20*y30<0,  y10=y30; r1=r3; end 
  end
  
  if abs(y10)<abs(y20) 
      r(i)=r1; 
  else 
      r(i)=r2; 
  end 

end


%---------------%
%  GRAPHICS     %
%---------------%
rccn=(rccn(2:end)+rccn(1:end-1))/2;  %middle of the size bins in um
dr=.5*(r(2:end)-r(1:end-1))*1e6;
r0=.5*(r(2:end)+r(1:end-1));  % output vector on droplet radii
r=r0;   %output droplet radii 

if k==1
    titl='marine CCN';
elseif k==2
    titl='clean continental CCN';
elseif k==3
    titl='background CCN';
elseif k==4
    titl='urban CCN';
end

figure
loglog(1e6*rccn,nccn);     hold on 
loglog(1e6*r0,nccn.*drccn./dr,'c') 
xlabel('dry CCN radius/droplet radius (um)') 
ylabel('CCN size distribution (cm^{-3} um^{-1})') 
legend('dry CCN','wetted CCN')
title(titl)

figure
loglog(1e6*rccn,1e6*r0) 
xlabel('dry CCN diameter (um)') 
ylabel('Equilibrium radius diameter (um)') 
title(['S=',num2str(S*100),'%'])

end

