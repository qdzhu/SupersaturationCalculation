classdef explore_GOAMAZON
    properties(Constant = true)
    end
    methods(Static)
        function SS = cal_ss(temp,alt,nconc,rmean,vel,opt)
            SS = zeros(size(rmean));
            [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=Constant; 
           for i = 1:numel(rmean)
               if nconc(i) == 0
                   continue;
               end
%                 aver_r = nansum(nccn.*drccn.*rccn_mn)/nansum(nccn.*drccn); % average radius
                %%% calculate temp, p, latent heat
                Tc = temp(i);
                T = To+Tc; % temperature in Kelvin
                Ho = alt(i);
                P=Po*exp(-g*Ho/(Ra*T)); % initial pressure
                aver_r = rmean(i);
                switch opt
                    case 'sim-var-dl'
                        L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
                        D=(2.26e-5+1.5e-7*Tc)*Po/P;  % diffusion coeff.
                    case 'sim-fixed-dl'
                        L = 2501000;
                        D = 0.23e-4;
                    case 'detailed'
                        L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
                        D=(2.26e-5+1.5e-7*Tc)*Po/P;  % diffusion coeff.
                        pa = P/(Ra*T); %%%% densityr for dry air;
                        Ew = (a0+Tc*(a1+Tc*(a2+Tc*(a3+Tc*(a4+Tc*(a5+Tc*a6))))))*100;  %vapor pressure over flat water surface 
                        K = 2.424e-2+7.95e-5*Tc;       % thermal conductivity
                        A0 = (g*L/(Cpa*Rv*T^2)-g/(Ra*T));
                        A6 = L^2/(Cpa*Rv*T^2);
                        A4 = P*Rv/(Ew*Ra);
                        Aw = (pl*L^2/(K*Rv*T^2)+pl*Rv*T/(Ew*D))^(-1);
                        Bw = 4*pi*pl*Aw/pa;
                        A = -A6*Bw*nconc(i)*aver_r;
                        B = (A4+A6)*Bw*nconc(i)*aver_r-A0*vel(i);
                        C = A0*vel(i);
                        SS(i) = (B-sqrt(B^2-4*A*C))/(2*A)*100;
                        continue;
                end
                
                A = (g*L/(Cpa*Rv*T^2)-g/(Ra*T));
                SS(i) = A*vel(i)/(4*pi*D*aver_r*nconc(i))*100;% quasi-state supersaturation
                %%A1=(9.806/(T+273.16)) * ( (2501000./(1005.6*461.53*(T+273.16) ) ) -(1/287.05) )
                %%SS=100.*A1*Wvel/((4.*3.14*0.23e-4*CDNC*rmean)) ;ss
            end
        end
        function read_goamazon_data()
            met_filepath = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/mei-iwg1';
            met_filepattern = 'aaf.iwg1001s.g1.goamazon.*';
            met_filedir = dir(fullfile(met_filepath,met_filepattern));
            dsd_filepath = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/comstock-fcdp';
            dsd_filepattern = 'FCDP_G1_*';
            dsd_filedir = dir(fullfile(dsd_filepath,dsd_filepattern));
%             r_range = [2.5:1:13.5,15:2:49];
            r_range = [1.5/2:1.5:(4.5+6)/2,7:2:17,(18+21)/2:3:(39+42)/2,44,48,100];
            for i=2:numel(met_filedir)
                metfile = fullfile(met_filepath,met_filedir(i).name);
                fileid = fopen(metfile);
                strfmt = strcat('%s %s',repmat(' %f ',1,44));
                metdata = textscan(fileid,strfmt,'Delimiter',',','headerLines', 2);
                fclose(fileid);
                date = strcat(extractBetween(met_filedir(i).name,26,29),'-',...
                   extractBetween(met_filedir(i).name,30,31),'-',...
                   extractBetween(met_filedir(i).name,32,33));
                time = datenum(string(metdata{2}));
                juliantime = (time-datenum(string(date)))*3600*24;
                lat = double(metdata{3});
                lon = double(metdata{4}); 
                alt = double(metdata{5});%m
                vel = double(metdata{29});
                vel(vel==-9999) = nan;
                incloudflag = double(metdata{37});
                phaseflag = double(metdata{38});
                temp = double(metdata{21});
                pres = double(metdata{24}); 
                
                dsdfile = fullfile(dsd_filepath,dsd_filedir(i).name);
                dsd = csvread(dsdfile,66,0);
                nconc = dsd(:,2)*1e-3;
                if size(dsd,1)==numel(metdata{1})
                    disp('dsd and met are one-to-one match');
                    n = size(dsd,1);
                else
                    n = min(numel(nconc),numel(lat));
                    nconc = nconc(1:n);
                    time = time(1:n);
                    juliantime = time(1:n);
                    vel = vel(1:n);
                    incloudflag = incloudflag(1:n);
                    temp = temp(1:n);
                    alt = alt(1:n);
                    pres = pres(1:n);
                end
                
                
                rmean = zeros(1,n);
                for j=1:min(numel(nconc),numel(lat))
                    if nconc(j)>0
                        rmean(j) = sum(dsd(j,3:end).*r_range)/nconc(j);
                    end
                    
                end
                opt = 'sim-fixed-dl';
                ss = explore_GOAMAZON.cal_ss(temp,alt,nconc,rmean,vel,opt);
                indx = incloudflag >0;
                figure;
                subplot(1,2,1);
                scatter(ss,vel,[],nconc,'filled');
                 h=colorbar;
                xmin = -5;
                xmax = 5;
                ymin = -15;
                ymax = 15;
                xlim([xmin,xmax]);
                ylim([ymin,ymax]);
                xlabel('SS,%');
                ylabel('W, ms^{-1}');
                ylabel(h,'Nc, cm^{-3}');
                caxis([0,40]);
                line([0,0],[ymin,ymax],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
                line([xmin,xmax],[0,0],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
                title(date);
                set(gca,'FontSize', 14);
                
                subplot(1,2,2);
                scatter(ss(indx),vel(indx),[],nconc(indx),'filled');
                 h=colorbar;
                xmin = -5;
                xmax = 5;
                ymin = -15;
                ymax = 15;
                xlim([xmin,xmax]);
                ylim([ymin,ymax]);
                xlabel('SS,%');
                ylabel('W, ms^{-1}');
                ylabel(h,'Nc, cm^{-3}');
                caxis([0,40]);
                line([0,0],[ymin,ymax],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
                line([xmin,xmax],[0,0],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
                title(date);
                set(gca,'FontSize', 14);
            end
            
        end
        function read_aerosol()
            %%% read midpoint of diameter bins
            dia_file = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/jimenez-smps/diameter-midpoint';
            fileid = fopen(dia_file);
            strfmt = repmat('%f ',1,90);
            dia = textscan(fileid,strfmt, 'Delimiter',' ');
            dia = cell2mat(dia);
            d_dia = dia(8:end)-dia(7:end-1);
            m_dia = dia(8:end);
            fclose(fileid);
            
            aerosol_file = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/jimenez-smps/CU GoA IOP1 SMPS dNdlogDm 26 Feb 2016 Level 2.txt';
            fileid = fopen(aerosol_file);
            strfmt = repmat('%f ',1,85);
            aerdata = textscan(fileid,strfmt,'Delimiter',' ','headerLines', 14);
            aerdata = cell2mat(aerdata);
            nm = aerdata(:,3:end);
            sttime = datenum('Jan 1, 1904')+aerdata(:,1)/86400-4/24; %% utc converts to local time
            edtime = datenum('Jan 1, 1904')+aerdata(:,2)/86400-4/24; %% utc converts to local time
            fclose(fileid);
            
            aerosol_file = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/jimenez-smps/CU GoA IOP1 and IOP2 SMPS integrated Number Surface Area Volume 26 Feb 2016 Level 2.txt';
            fileid = fopen(aerosol_file);
            strfmt = repmat('%f ',1,5);
            aerdata = textscan(fileid,strfmt,'Delimiter',' ','headerLines', 14);
            aerdata = cell2mat(aerdata);
            nm = aerdata(:,3);
            aer_sttime = datenum('Jan 1, 1904')+aerdata(:,1)/86400-4/24; %% utc converts to local time
            aer_edtime = datenum('Jan 1, 1904')+aerdata(:,2)/86400-4/24; %% utc converts to local time
            fclose(fileid);
            
        end
        function reanalysis_goamazon_data()
            %%% read aerosol data from Fan's paper
            fileid = fopen('/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/jimenez-smps/aerosol-number-conc'); 
            strfmt = '%s %f %s %f %f %f %f';
            aerdata = textscan(fileid,strfmt,'Delimiter',',','headerLines', 1);
            aerosol_date = datenum(aerdata{1});
            aerosol_nm = aerdata{2};
            aerosol_flag = aerdata{5};
            start_conv = aerdata{6};
            end_conv = aerdata{7};
            
            %%% read rwp data and reflectivity map
            rwp_filepath = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/giangrande-rwp';
            sipam_filepath = '/Users/monicazhu/Box/AtmosDynamic/GOAMAZON/schumacher-sband_radar';
            sipam_filepattern = '2014*';
            sipam_filedir = dir(fullfile(sipam_filepath,sipam_filepattern));
            
            tot_aer = [];
            tot_hgt = [];
            tot_vel = [];
            tot_ref = [];
            tot_flg = [];
            tot_sipam_radar = [];
            tot_sipam_vel = [];
            sipam_height = [];
            sipam_lon = [];
            sipam_lat = [];
            for i=1:numel(sipam_filedir)
                thisname = sipam_filedir(i).name;
                thisdate = strcat(extractBetween(thisname,1,4),'-',extractBetween(thisname,5,6),'-',extractBetween(thisname,7,8));
                thisaerosol = aerosol_nm(aerosol_date == datenum(thisdate));
                thisflg = aerosol_flag(aerosol_date == datenum(thisdate));
                %%% read rwp data
                rwppattern = strcat('maorwpcls.',thisname,'.cdf');
                rwpfiledir = dir(fullfile(rwp_filepath,rwppattern));
                rwpfile = ncinfo(fullfile(rwpfiledir.folder,rwpfiledir.name));
                rwp_time = ncread(rwpfile.Filename,'time_offset');
                rwp_height = ncread(rwpfile.Filename,'height');
                rwp_ref = ncread(rwpfile.Filename,'ReflectivityUAZR');
                rwp_w = ncread(rwpfile.Filename,'VerticalVelocity');
                rwp_flag = ncread(rwpfile.Filename,'EchoClassification');
                %%% local time 1100-1900
                rwp_ref = rwp_ref(:,15*600+1:23*600);
                rwp_ref(rwp_ref<-100) = nan;
                rwp_w = rwp_w(:,15*600+1:23*600);
                rwp_w(rwp_w<-100) = nan;
                rwp_flag = rwp_flag(:,15*600+1:23*600);
                rwp_time = rwp_time(15*600+1:23*600);
                
%                 conv_filter= find(rwp_flag(10,:) ==6);
%                 conv_time = conv_filter/600+15;
                %%% find out convection events
%                 rwp_ref(rwp_flag~=6) = nan;
%                 rwp_w(rwp_flag~=6) = nan;
                ext_w = zeros(size(rwp_height));
                ext_ref = zeros(size(rwp_height));
                indx = ~(any(rwp_flag(83:end,:))>0);
                rwp_ref(:,indx) = nan;
                rwp_w(:,indx) = nan;
                conv_time = find(indx == 0)/600+15;
%                 
                %determine when the convection is by eye
%                 indx = rwp_time/3600 > start_conv(i)+4 & rwp_time/3600<end_conv(i)+4;
%                 rwp_ref(:,~indx) = nan;
%                 rwp_w(:,~indx) = nan;
                
%                 disp(thisdate);
%                 disp(min(conv_time));
%                 disp(max(conv_time));
%                 
                for j=1:numel(rwp_height)
                    this_w = rwp_w(j,:);
                    this_ref = rwp_ref(j,:);
                    this_ref = this_ref(this_w>0);
                    this_w = this_w(this_w>0);
                    indx = this_w > prctile(this_w,90);
                    if sum(indx)>5
                        ext_w(j) = nanmean(this_w(indx));
                        ext_ref(j) = nanmean(this_ref(indx));
                    else
                        ext_w(j) = nan;
                        ext_ref(j) = nan;
                    end
                end
                
                tot_aer = [tot_aer,zeros(size(rwp_height))+thisaerosol];
                tot_flg = [tot_flg,zeros(size(rwp_height))+thisflg];
                tot_hgt = [tot_hgt,rwp_height];
                tot_ref = [tot_ref,ext_ref];
                tot_vel = [tot_vel,ext_w];
                
%                 figure;
%                 [meshtime,meshhgt] = meshgrid(rwp_time(15*600+1:23*600)/3600-4,rwp_height);
%                 s = subplot(2,1,1);
%                 pcolor(meshtime,meshhgt,rwp_w);
%                 colormap(s,'blue_red_cmap');
%                 shading flat;
%                 h = colorbar;
%                 ylabel(h,'Vertical Velocity (m/s)');
%                 caxis([-15,15]);
%                 xlabel('hour of day (local)');
%                 ylabel('Height (km)');
%                 title(thisname);
%                 set(gca,'FontSize', 14);   
%                 
%                 s = subplot(2,1,2);
%                 pcolor(meshtime,meshhgt,rwp_ref);
%                 colormap(s,'jet');
%                 shading flat;
%                 h = colorbar;
%                 ylabel(h,'Radar reflectivity (dBZ)');
%                 caxis([0,50]);
%                 xlabel('hour of day (local)');
%                 ylabel('Height (km)');
%                 set(gca,'FontSize', 14);   

                %%% read sipam
                sipamdir = dir(fullfile(sipam_filepath,thisname,strcat('*',thisname,'*')));
                sipamtime = conv_time_file(sipamdir);
               
                sipam_radar = [];
                
                sipam_vel = [];
                exist_flg = 0;
                for k=1:numel(sipamdir)
%                     if (sipamtime(k)-datenum(thisdate))*24 > start_conv(i)+4 &&  (sipamtime(k)-datenum(thisdate))*24<end_conv(i)+4;
                    if min(abs(sipamtime(k)-datenum(thisdate)-conv_time/24))< 12/(60*24) %&& sipamtime(k) < datenum(thisdate)+23/24
%                     if sipamtime(k) > datenum(thisdate)+15/24 && sipamtime(k) < datenum(thisdate)+23/24
                        sipam_info = ncinfo(fullfile(sipamdir(k).folder,sipamdir(k).name));
                        if isempty(sipam_lon)
                            sipam_lon = ncread(sipam_info.Filename,'lon0');
                            sipam_lat = ncread(sipam_info.Filename,'lat0');
                            sipam_height = ncread(sipam_info.Filename,'z0');
                        end
                        this_radar = ncread(sipam_info.Filename,'DBZc');
                        this_radar(this_radar<-100) = nan;
                        this_vel = ncread(sipam_info.Filename,'VEL');
                        this_vel(this_vel<-100) = nan;
                        sipam_radar = cat(4,sipam_radar,this_radar);
                        sipam_vel = cat(4,sipam_vel,this_vel);
                        exist_flg = 1;
                    end
                    if k==numel(sipamdir) && exist_flg ==0
                        sipam_radar = nan(size(tot_sipam_radar,1),size(tot_sipam_radar,2),size(tot_sipam_radar,3),1);
                        sipam_vel = nan(size(tot_sipam_radar,1),size(tot_sipam_radar,2),size(tot_sipam_radar,3),1);
                    end
                end
                sipam_radar(sipam_radar<prctile(sipam_radar,90,4)) = nan;
                sipam_vel(sipam_vel<prctile(sipam_vel,90,4)) = nan;
                ext_sipam_radar = nanmean(sipam_radar,4);
                ext_sipam_vel = nanmean(sipam_vel,4);
                tot_sipam_radar = cat(4,tot_sipam_radar,ext_sipam_radar);
                tot_sipam_vel = cat(4,tot_sipam_vel,ext_sipam_vel);
            end
            save('goamazon_90th_conv_fan','tot_aer','tot_hgt','tot_ref','tot_vel','tot_flg','sipam_lon','sipam_lat','sipam_height','tot_sipam_radar','tot_sipam_vel');

        end
        function plot_goamazon_data()
            data = load('/Users/monicazhu/Box/AtmosDynamic/goamazon_90th_conv_fan.mat');
            aer = data.tot_aer;
            hgt = data.tot_hgt;
            vel = data.tot_vel;
            ref = data.tot_ref;
            flg = data.tot_flg;
            sipam_radar = data.tot_sipam_radar;
            sipam_vel = data.tot_sipam_vel;
            sipam_height = data.sipam_height;
            sipam_vel = data.tot_sipam_vel;
            lon = data.sipam_lon;
            lat = data.sipam_lat;
            %%% plot sipam data
            sipam_radar(sipam_radar<0) = nan;
            indx = find(flg(1,:) ==0);
            clean_radar = nanmean(sipam_radar(:,:,:,indx),4);
            indx = find(flg(1,:) ==1);
            sc_radar = nanmean(sipam_radar(:,:,:,indx),4);
            indx = find(flg(1,:) ==2);
            su_radar = nanmean(sipam_radar(:,:,:,indx),4);
            indx = find(flg(1,:) ==4);
            urban_radar = nanmean(sipam_radar(:,:,:,indx),4);
            
            figure;
%             pos3 = [-60.6,-3.2];
            pos3 = [-60.55,-3.19];
            [x,y] = find(lon >= pos3(1)-0.01 & lon<= pos3(1)+0.01 &...
                lat >= pos3(2)-0.01 & lat <= pos3(2)+0.01);
            rwp_clean = squeeze(nanmean(nanmean(clean_radar(x,y,:),1),2));
            rwp_sc = squeeze(nanmean(nanmean(sc_radar(x,y,:),1),2));
            rwp_su = squeeze(nanmean(nanmean(su_radar(x,y,:),1),2));
            rwp_urban = squeeze(nanmean(nanmean(urban_radar(x,y,:),1),2));
            
            l(1) = line(rwp_clean,sipam_height,'linestyle','-','color','k');
            l(2) = line(rwp_sc,sipam_height,'linestyle','-','color','g');
            l(3) = line(rwp_su,sipam_height,'linestyle','-','color','b');
            l(4) = line(rwp_urban,sipam_height,'linestyle','-','color','r');
            legend(l,[{'500-1000'},{'1000-1900'},{'1900-3000'},{'>3000'}]);
            xlabel('Reflectivity (dBZ)');
            ylabel('Height (km)');
            set(gca,'FontSize', 14); 
            
            figure;
            subplot(2,2,1);
            pcolor(lon,lat,squeeze(clean_radar(:,:,16)));
            shading interp;
            h = colorbar;
%             h = colorbar;
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
%             colormap('jet');
            caxis([0,40]);
            title('500-1000');
            line(-60,-3.1,'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            line(-60.6,-3.2,'linestyle','none','marker','o','color','b','linewidth',2,'MarkerSize',5);
            line(-60.55,-3.19,'linestyle','none','marker','o','color','r','linewidth',2,'MarkerSize',5);
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            set(gca,'FontSize', 14); 
            
            subplot(2,2,2);
            pcolor(lon,lat,squeeze(sc_radar(:,:,16)));
            shading interp;
            h = colorbar;
            caxis([0,40]);
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
%             colormap('jet');
            title('1000-1900');
            line(-60,-3.1,'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            line(-60.6,-3.2,'linestyle','none','marker','o','color','b','linewidth',2,'MarkerSize',5);
            line(-60.55,-3.19,'linestyle','none','marker','o','color','r','linewidth',2,'MarkerSize',5);
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            set(gca,'FontSize', 14); 
            
            subplot(2,2,3);
            pcolor(lon,lat,squeeze(su_radar(:,:,16)));
            shading interp;
            h = colorbar;
            caxis([0,40]);
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
%             colormap('jet');
            title('1900-3000');
            line(-60,-3.1,'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            line(-60.6,-3.2,'linestyle','none','marker','o','color','b','linewidth',2,'MarkerSize',5);
            line(-60.55,-3.19,'linestyle','none','marker','o','color','r','linewidth',2,'MarkerSize',5);
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            set(gca,'FontSize', 14); 
            
            subplot(2,2,4);
            pcolor(lon,lat,squeeze(urban_radar(:,:,16)));
            shading interp;
            h = colorbar;
            caxis([0,40]);
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
            line(-60,-3.1,'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            line(-60.6,-3.2,'linestyle','none','marker','o','color','b','linewidth',2,'MarkerSize',5);
            line(-60.55,-3.19,'linestyle','none','marker','o','color','r','linewidth',2,'MarkerSize',5);
%             colormap('jet');
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            title('>3000');
            set(gca,'FontSize', 14); 
            
            figure;
            s=subplot(1,3,1);
            pcolor(lon,lat,squeeze(sc_radar(:,:,16))-squeeze(clean_radar(:,:,16)));
            shading interp;
            h = colorbar;
%             h = colorbar;
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
            colormap(s,blue_red_cmap);
%             colormap('jet');
            caxis([-15,15]);
            title('500-1000 vs 1000-1900');
            line([-60.6,-60],[-3.2,-3.1],'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            set(gca,'FontSize', 14); 
            
            s=subplot(1,3,2);
            pcolor(lon,lat,squeeze(su_radar(:,:,16))-squeeze(sc_radar(:,:,16)));
            shading interp;
            h = colorbar;
            caxis([-15,15]);
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
            colormap(s,blue_red_cmap);
%             colormap('jet');
            title('1900-3000 vs 1000-1900');
            line([-60.6,-60],[-3.2,-3.1],'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            set(gca,'FontSize', 14); 
            
            s = subplot(1,3,3);
            pcolor(lon,lat,squeeze(urban_radar(:,:,16))-squeeze(su_radar(:,:,16)));
            shading interp;
            h = colorbar;
            caxis([-15,15]);
            ylabel(h,'Reflectivity (dBZ)');
            xlabel('longitude');
            ylabel('Latitude');
            colormap(s,blue_red_cmap);
            title('>3000 vs 1900-3000');
            line([-60.6,-60],[-3.2,-3.1],'linestyle','none','marker','d','color','k','linewidth',2,'MarkerSize',12);
            xlim([-61,-60]);
            ylim([-3.5,-2.9]);
            set(gca,'FontSize', 14); 
            %%% plot aerosol concentration - height - vertical velocity
            figure;
            subplot(2,1,1);
            st_aer = reshape(aer,17*140,1);
            st_hgt = reshape(hgt,17*140,1);
            st_vel = reshape(vel,17*140,1);
            scatter(st_aer,st_hgt,200,st_vel,'marker','.');
            h = colorbar;
            caxis([0,15]);
            ylim([0,16]);
            ylabel(h,'Vertical velocity (m/s)');
            ylabel('Height (km)');
            xlabel('Aerosol concentration (D>15nm; cm^{-3})');
            set(gca,'FontSize', 14); 
            
            subplot(2,1,2);
            st_ref = reshape(ref,17*140,1);
            scatter(st_aer,st_hgt,200,st_ref,'marker','.');
            h = colorbar;
            ylabel(h,'Reflectivity (dBZ)');
            caxis([0,50]);
            ylim([0,16]);
            ylabel('Height (km)');
            xlabel('Aerosol concentration (D>15nm; cm^{-3})');
            set(gca,'FontSize', 14); 
            
            cleanindx = find(flg(1,:) ==0);
            clean_w = nanmean(vel(:,cleanindx),2);
            [clean_w, clean_std] = bin_omisp_height(clean_w,squeeze(hgt(:,1)),sipam_height);
            clean_ref = nanmean(ref(:,cleanindx),2);
            [clean_ref, clean_ref_std] = bin_omisp_height(clean_ref,squeeze(hgt(:,1)),sipam_height);
            
            scindx = find((flg(1,:)) == 1);
            sc_w = nanmean(vel(:,scindx),2);
            [sc_w, sc_std] = bin_omisp_height(sc_w,squeeze(hgt(:,1)),sipam_height);
            sc_ref = nanmean(ref(:,scindx),2);
            [sc_ref, sc_ref_std] = bin_omisp_height(sc_ref,squeeze(hgt(:,1)),sipam_height);
            
            suindx = find(flg(1,:) == 2);
            su_w = nanmean(vel(:,suindx),2);
            [su_w, su_std] = bin_omisp_height(su_w,squeeze(hgt(:,1)),sipam_height);
            su_ref = nanmean(ref(:,suindx),2);
            [su_ref, su_ref_std] = bin_omisp_height(su_ref,squeeze(hgt(:,1)),sipam_height);
            
            urbanindx = find(flg(1,:) ==4);
            urban_w = nanmean(vel(:,urbanindx),2);
            [urban_w, urban_std] = bin_omisp_height(urban_w,squeeze(hgt(:,1)),sipam_height);
            urban_ref = nanmean(ref(:,urbanindx),2);
            [urban_ref, urban_ref_std] = bin_omisp_height(urban_ref,squeeze(hgt(:,1)),sipam_height);
            
            figure;
            subplot(1,2,1);
            hold on;
            fill_x = [clean_w-clean_std,flip(clean_w+clean_std)];
            fill_h = cat(1,sipam_height,flip(sipam_height));
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'k','Linestyle','none');
            alpha(f,.1);
            
            fill_x = [sc_w-sc_std,flip(sc_w+sc_std)];
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'g','Linestyle','none');
            alpha(f,.1);
            
            fill_x = [su_w-su_std,flip(su_w+su_std)];
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'b','Linestyle','none');
            alpha(f,.1);
            
            fill_x = [urban_w-urban_std,flip(urban_w+urban_std)];
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'r','Linestyle','none');
            alpha(f,.1);
            
            l(1) = line(clean_w,sipam_height,'linestyle','-','color','k');
            l(2) = line(sc_w,sipam_height,'linestyle','-','color','g');
            l(3) = line(su_w,sipam_height,'linestyle','-','color','b');
            l(4) = line(urban_w,sipam_height,'linestyle','-','color','r');
            legend(l,[{'500-1000'},{'1000-1900'},{'1900-3000'},{'>3000'}]);
            xlabel('Veltical velocity (m/s)');
            ylabel('Height (km)');
            
            subplot(1,2,2);
            hold on;
            fill_x = [clean_ref-clean_ref_std,flip(clean_ref+clean_ref_std)];
            fill_h = cat(1,sipam_height,flip(sipam_height));
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'k','Linestyle','none');
            alpha(f,.1);
            
            fill_x = [sc_ref-sc_ref_std,flip(sc_ref+sc_ref_std)];
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'g','Linestyle','none');
            alpha(f,.1);
            
            fill_x = [su_ref-su_ref_std,flip(su_ref+su_ref_std)];
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'b','Linestyle','none');
            alpha(f,.1);
            
            fill_x = [urban_ref-urban_ref_std,flip(urban_ref+urban_ref_std)];
            fill_y = fill_h(~isnan(fill_x));
            fill_x = fill_x(~isnan(fill_x));
            f = fill(fill_x,fill_y, 'r','Linestyle','none');
            alpha(f,.1);
            
            l(1) = line(clean_ref,sipam_height,'linestyle','-','color','k');
            l(2) = line(sc_ref,sipam_height,'linestyle','-','color','g');
            l(3) = line(su_ref,sipam_height,'linestyle','-','color','b');
            l(4) = line(urban_ref,sipam_height,'linestyle','-','color','r');
            legend(l,[{'500-1000'},{'1000-1900'},{'1900-3000'},{'>3000'}]);
            xlabel('Reflectivity (dBZ)');
            ylabel('Height (km)');
        end
        
            
    end
end

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
g= 9.81;        % acceleration of gravity
k_mu=sqrt(1/(2*pi)); 
k_ml=4*pi*pl/3; 
end

function sipamtime = conv_time_file(sipam)
    n = numel(sipam);
    sipamtime = zeros(n,1);
    for i=1:n
        thisname = strsplit(sipam(i).name,'_');
        thisdate = strcat(extractBetween(thisname(3),1,4),'-',extractBetween(thisname(3),5,6),'-',extractBetween(thisname(3),7,8));
        thistime = strcat(extractBetween(thisname(4),1,2),':',extractBetween(thisname(4),3,4),':',extractBetween(thisname(4),5,6));
        sipamtime(i) = datenum(strcat(thisdate,{' '},thistime));          
    end
end

function [bin_values,bin_errors] = bin_omisp_height(data_vals,height,bin_midpoints)
    binmode = 'mean';
    delta = diff(bin_midpoints);
    bins = zeros(numel(bin_midpoints,3));
    for a=1:numel(bin_midpoints)
        % The first and last bins need special handling, since diff produces a
        % vector one shorter than its input
        if a==1
            D = delta(1);
        elseif a==numel(bin_midpoints)
            D = delta(end);
        else
            D = max(delta((a-1):a));
        end
        % bins will have three columns: the bottom of the bin, the bin
        % midpoint, and the top of the top.
        bins(a,1) = bin_midpoints(a) - D/2;
        bins(a,2) = bin_midpoints(a);
        bins(a,3) = bin_midpoints(a) + D/2;
    end
    bin_values = zeros(1,size(bins,1));
    bin_errors = zeros(1,size(bins,1));
    for a=1:numel(bin_values)
        bin_data_vals = data_vals(height > bins(a,1) & height <= bins(a,3));
        if strcmpi(binmode,'mean')
            bin_values(a) = nanmean(bin_data_vals(:));
            bin_errors(a) = nanstd(bin_data_vals(:)) / sqrt(numel(bin_data_vals));
        elseif strcmpi(binmode, 'median')
            bin_values(a) = nanmedian(bin_data_vals(:));
            bin_errors(:,a) = quantile(bin_data_vals(:),[0.1,0.9]);
    %         bin_errors(:,a) = quantile(bin_data_vals(:),[0.25,0.75]);
        else
            bin_values{a} = bin_data_vals(:);
        end
    end 

end