classdef recalculate_CAIPEEX_result
    properties(Constant = true)
    end
    methods(Static)
        % a single case is used to get familar with the CAIPEEX data,
        % especially the cloud size distribution data without any
        % documentation for units.
        function single_case()
            [Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,...
            Mms,alpha,w,Po,To,g]=Constant;  
            T = To-10.87; % temperature in Kelvin
            Tc = -10.87; % temperature in Celsius degress from CAIPEEX dataset
            Ho = 2904;  % altitude from Figure
            P=Po*exp(-g*Ho/(Ra*T)); % initial pressure
            % CONSTANTS TEMPERATURE & PRESSURE DEPENDENT 
            D=(2.26e-5+1.5e-7*Tc)*Po/P;  % diffusion coeff.
            L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
            K=2.424e-2+7.95e-5*Tc;       % thermal conductivity
            sigma=7.564e-2-1.43e-4*Tc;   % water surface tension

            filename = 'Taraha_ex_PSD.csv';
            data = csvread(filename,6,0);
            rccn = data(:,1);
            sp = 10.^data(:,2); %DSD distribution from plot digitizer
            drccn = rccn(2:end)-rccn(1:end-1); %  %length of the size bins
            nccn = (sp(1:end-1)+sp(2:end))/2; %middle of the size bins
            rccn_mn = (rccn(1:end-1)+rccn(2:end))/2;% mean of size
            %%%Nconc = 38 cm^-3
            Nconc=nansum(nccn.*drccn); % droplet number concentration 
            aver_r = nansum(nccn.*drccn.*rccn_mn)/nansum(nccn.*drccn); % average radius
            A = (g*L/(Cpa*Rv*T^2)-g/(Ra*T));
            Vel = 3.36;% vertical velocity
            S = A*Vel/(4*pi*D*aver_r*Nconc);% quasi-state supersaturation

        end
        function plot_dsd_single_case()
            filepath = '/Users/monicazhu/Documents/MATLAB/AtmosDynamic';
            metfile = '200906220753_1Hz.csv';
            dsdfile = '20090622_30000_1Hz_dsd.csv';
            dsddata = csvread(fullfile(filepath,dsdfile),9,1);
            indx = dsddata(:,1) < 50;
            dsddata = dsddata(indx,:);
            diameter = dsddata(:,1);
            dsd = dsddata(:,2:end);
            
            metdata = csvread(fullfile(filepath,metfile),1,0);
            juliantime = metdata(:,1);
            indx = juliantime == 95815;
            example_dsd = dsd(:,indx);
            example_dsd = example_dsd./((3.14/6.)*(diameter*1e-4).^3)*1e-6;
            
            filename = 'Taraha_ex_PSD.csv';
            data = csvread(filename,6,0);
            rccn = data(:,1);
            sp = 10.^data(:,2); %DSD distribution from plot digitizer
            
            figure;
            subplot(1,2,1);
            hold on;
            plot(diameter,example_dsd,rccn,sp);
            legend('From Dataset','From paper');
            xlabel('D ({\mu}m)')
            ylabel('dN/dD (cm^{-3}{\mu}m^{-1})');
            set(gca, 'YScale', 'log')
            hold off;
            
%             subplot(1,2,2);
%             plot(diameter,rescale_sp./example_dsd);
%             ylabel('Paper/Dataset');
%             xlabel('D ({\mu}m)')
            
            
        end
        % the utility function reading inputs from a daily file and
        % calculating supersaturation.
        function [SS,nconc,vel,alt,sec,lat,lon,lwc,af,rm]=Daily_cal(metfile,dsdfile,opt)
            filepath = '/Users/monicazhu/Box/AtmosDynamic/CAIPEEX';
            %%% read julian time, altitude, temperature, vertical velocity from met datafile
            metdata = csvread(fullfile(filepath,metfile),1,0);
            juliantime = metdata(:,1);
            sec = metdata(:,2);
            temp = metdata(:,3);%%% celcius degree
            lat = metdata(:,8);
            lon = metdata(:,9);
            alt = metdata(:,10);%%% m
            vel = metdata(:,15);%%% m/s
            nconc = metdata(:,20);%%% #cm^-3
            lwc = metdata(:,23);
            af = metdata(:,24);
            dsddata = csvread(fullfile(filepath,dsdfile),9,1);
            indx = dsddata(:,1) <= 50;
            dsddata = dsddata(indx,:);
            diameter = dsddata(:,1);
            dsd = dsddata(:,2:end);
            dsd = dsd./((3.14/6.)*(diameter*1e-4).^3)*1e-6;
            
            if size(dsd,2) == numel(juliantime)
                disp('dsd data is 1-to-1 match to met data');
            end
            [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=recalculate_CAIPEEX_result.Constant; 
            SS = zeros(size(juliantime));
            rm = zeros(size(juliantime));
            for i = 1:numel(juliantime)
                sp = dsd(:,i);
                drccn = diameter(2:end)-diameter(1:end-1);
                nccn = (sp(1:end-1)+sp(2:end))/2; %middle of the size bins
                rccn_mn = (diameter(1:end-1)+diameter(2:end))/2;% mean of size
                indx = rccn_mn <2000;
                aver_r = nansum(nccn(indx).*drccn(indx).*rccn_mn(indx))/nansum(nccn(indx).*drccn(indx)); % average radius
                rm(i) = aver_r;
                %%% calculate temp, p, latent heat
                Tc = temp(i);
                T = To+Tc; % temperature in Kelvin
                Ho = alt(i);
                P=Po*exp(-g*Ho/(Ra*T)); % initial pressure
%                 
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
        % plot procedure
        % explore the relationship between vertical wind velocity
        function plot_vel_ss()
            opt = 'sim-fixed-dl';
            metfile = '200906160818_1Hz.csv';
            dsdfile = '20090616_30000_1Hz_dsd.csv';
            [SS.d1,nconc.d1,vel.d1,alt.d1,sec.d1,lat.d1,lon.d1,lwc.d1,af.d1]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200906210758_1Hz.csv';
            dsdfile = '20090621_30000_1Hz_dsd.csv';
            [SS.d2,nconc.d2,vel.d2,alt.d2,sec.d2,lat.d2,lon.d2,lwc.d2,af.d2]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200906220753_1Hz.csv';
            dsdfile = '20090622_30000_1Hz_dsd.csv';
            [SS.d3,nconc.d3,vel.d3,alt.d3,sec.d3,lat.d3,lon.d3,lwc.d3,af.d3]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908180730_1Hz.csv';
            dsdfile = '20090818_30000_1Hz_dsd.csv';
            [SS.d4,nconc.d4,vel.d4,alt.d4,sec.d4,lat.d4,lon.d4,lwc.d4,af.d4]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908230850_1Hz.csv';
            dsdfile = '20090823_30000_1Hz_dsd.csv';
            [SS.d5,nconc.d5,vel.d5,alt.d5,sec.d5,lat.d5,lon.d5,lwc.d5,af.d5]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908240858_1Hz.csv';
            dsdfile = '20090824_30000_1Hz_dsd.csv';
            [SS.d6,nconc.d6,vel.d6,alt.d6,sec.d6,lat.d6,lon.d6,lwc.d6,af.d6]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908250853_1Hz.csv';
            dsdfile = '20090825_30000_1Hz_dsd.csv';
            [SS.d7,nconc.d7,vel.d7,alt.d7,sec.d7,lat.d7,lon.d7,lwc.d7,af.d7]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            fns = {'d1','d2','d3','d4','d5','d6','d7'};
            titlestr = {'16 Jun','21 Jun','22 Jun','18 Aug','23 Aug','24 Aug','25 Aug'};

            figure;
            for a =1:numel(fns)
                subplot(3,3,a);
%                 scatter(SS.(fns{a}),vel.(fns{a}),2.^(alt.(fns{a})/1000),nconc.(fns{a}),'filled');
                scatter(SS.(fns{a}),vel.(fns{a}),[],lwc.(fns{a}),'filled');
                h=colorbar;
                xmin = -4;
                xmax = 7;
                ymin = -15;
                ymax = 15;
                xlim([xmin,xmax]);
                ylim([ymin,ymax]);
                xlabel('SS,%');
                ylabel('W, ms^{-1}');
                ylabel(h,'Nc, cm^{-3}');
%                 caxis([0,400]);
                line([0,0],[ymin,ymax],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
                line([xmin,xmax],[0,0],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
                disp(max(SS.(fns{a})));
                title(titlestr{a});
                set(gca,'FontSize', 14);
            end
            figure;
            for a =1:numel(fns)-4
                this_ss = SS.(fns{a});
                this_lwc = lwc.(fns{a});
                
                subplot(3,2,2*a-1);
                scatter(this_ss(this_ss>0),this_lwc(this_ss>0));
                subplot(3,2,2*a);
                scatter(this_ss(this_ss<0),this_lwc(this_ss<0));
            end
%             hold on;
%             bubsizes = 2.^[3,5,7];
%             scatter([2.2,2.2,2.2],[-12,-10,-8],bubsizes,'blue','filled');
%             %text([2.3,2.3,2.3],[-12,-10,-8],{'3km','5km','7km'});
%             rectangle('Position',[2.0,-13,1.8,7],'LineWidth',2);
%             hold off;
        
%             legend(legendtry);
            
        end
        % explore the relationship between adiabatic fraction and
        % supersaturation
        function explore_af_ss()
            opt = 'sim-fixed-dl';
            metfile = '200906160818_1Hz.csv';
            dsdfile = '20090616_30000_1Hz_dsd.csv';
            [SS.d1,nconc.d1,vel.d1,alt.d1,sec.d1,lat.d1,lon.d1,lwc.d1,af.d1,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200906210758_1Hz.csv';
            dsdfile = '20090621_30000_1Hz_dsd.csv';
            [SS.d2,nconc.d2,vel.d2,alt.d2,sec.d2,lat.d2,lon.d2,lwc.d2,af.d2,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200906220753_1Hz.csv';
            dsdfile = '20090622_30000_1Hz_dsd.csv';
            [SS.d3,nconc.d3,vel.d3,alt.d3,sec.d3,lat.d3,lon.d3,lwc.d3,af.d3,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908180730_1Hz.csv';
            dsdfile = '20090818_30000_1Hz_dsd.csv';
            [SS.d4,nconc.d4,vel.d4,alt.d4,sec.d4,lat.d4,lon.d4,lwc.d4,af.d4,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908230850_1Hz.csv';
            dsdfile = '20090823_30000_1Hz_dsd.csv';
            [SS.d5,nconc.d5,vel.d5,alt.d5,sec.d5,lat.d5,lon.d5,lwc.d5,af.d5,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908240858_1Hz.csv';
            dsdfile = '20090824_30000_1Hz_dsd.csv';
            [SS.d6,nconc.d6,vel.d6,alt.d6,sec.d6,lat.d6,lon.d6,lwc.d6,af.d6,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908250853_1Hz.csv';
            dsdfile = '20090825_30000_1Hz_dsd.csv';
            [SS.d7,nconc.d7,vel.d7,alt.d7,sec.d7,lat.d7,lon.d7,lwc.d7,af.d7,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            ss_all = cat(1,SS.d1,SS.d2,SS.d3,SS.d4,SS.d5,SS.d6,SS.d7);
            af_all = cat(1,af.d1,af.d2,af.d3,af.d4,af.d5,af.d6,af.d7);
            vel_all = cat(1,vel.d1,vel.d2,vel.d3,vel.d4,vel.d5,vel.d6,vel.d7);
            nconc_all = cat(1,nconc.d1,nconc.d2,nconc.d3,nconc.d4,nconc.d5,nconc.d6,nconc.d7);
            alt_all = cat(1,alt.d1,alt.d2,alt.d3,alt.d4,alt.d5,alt.d6,alt.d7);
            lwc_all = cat(1,lwc.d1,lwc.d2,lwc.d3,lwc.d4,lwc.d5,lwc.d6,lwc.d7);
            
            indx = ss_all <20;
            ss_all = ss_all(indx);
            af_all = af_all(indx);
            vel_all = vel_all(indx);
            lwc_all = lwc_all(indx);
            nconc_all = nconc_all(indx);
            alt_all = alt_all(indx);
            
            ss_cloud_cntr_all = ss_all(lwc_all >= prctile(lwc_all(:),5));
            ss_cloud_edge_all = ss_all(lwc_all < prctile(lwc_all(:),5));
            
            figure;
            hold on;
            histogram(ss_cloud_cntr_all,30,'Normalization','probability');
            histogram(ss_cloud_edge_all,30,'Normalization','probability');
            legend('Non-edge','edge');
            xlabel('SS (%)');
            ylabel('Probability');
            hold off;
           
            % pdf 
            figure;
            subplot(1,2,1);
            histogram(af_all(ss_all>0),'BinWidth',0.02);
            xlim([0,1]);
            ylabel('Contribution Fraction');
            xlabel('Adiabatic Fraction');
            set(gca,'FontSize', 14);
            subplot(1,2,2);
            histogram(af_all(ss_all<0),'BinWidth',0.02);
            xlim([0,1]);
            ylabel('Contribution Fraction');
            xlabel('Adiabatic Fraction');
            set(gca,'FontSize', 14);
            % ss vs vel, separate updraft and downdraft
            figure;
            scatter(ss_all,vel_all,[],nconc_all,'filled');
            h = colorbar;
            xmin = -4;
            xmax = 7;
            ymin = -15;
            ymax = 15;
            xlim([xmin,xmax]);
            ylim([ymin,ymax]);
            xlabel('SS,%');
            ylabel('W, ms^{-1}');
            ylabel(h,'Nconc cm^{-3}');
%                 caxis([0,400]);
            line([0,0],[ymin,ymax],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');
            line([xmin,xmax],[0,0],'color','k','linestyle',':','linewidth',2,'HandleVisibility','off');

            set(gca,'FontSize', 14);
                
            
            % histogram of vel categorized by ss
            vel_up = zeros(8,4);
            for i=1:4
                if i==4
                   vel_up_bin = vel_all(ss_all>2); 
                end
                if i==3
                    vel_up_bin = vel_all(ss_all<2 & ss_all>1);
                end
                if i==2
                    vel_up_bin = vel_all(ss_all<1 & ss_all>0.5);
                end
                if i==1
                    vel_up_bin = vel_all(ss_all<0.5 & ss_all>0);
                end
                vel_bin = 0:2:16;
                for j=1:8
                    indx = vel_up_bin >vel_bin(j) & vel_up_bin < vel_bin(j+1);
                    vel_up(j,i) = sum(indx)/numel(vel_up_bin);
                end
            end
            
            figure;
            b=bar(1:2:15,vel_up,'group');
%             labelstr = {'0-2','2-4','4-6','6-8','>8'};
            legend('SS <0.5','SS 0.5-1','SS 1-2','SS >2');
%             set(gca,'yticklabel',labelstr); 
            xlabel('W (m/s)');
            ylabel('Fraction');
            
            vel_up1 = vel_all(ss_all>2);
            vel_up2 = vel_all(ss_all<2 & ss_all>1);
            vel_up3 = vel_all(ss_all<1 & ss_all>0.5);
            vel_up4 = vel_all(ss_all<0.5 & ss_all>0);
            figure;
            histogram(vel_up1);
            
            % vel vs af vs SS
            figure;
            af_g = [0.02,0.05,0.1,0.5,1];
            vel_up_g = [0,2,4,6,8,20];
            vel_down_g =[0,-2,-4,-6,-8,-20];% [-20,-8,-6,-4,-2,0];
            ss_vel_up = zeros(numel(vel_up_g)-1,4);
            af_vel_up = zeros(numel(vel_up_g)-1,4);
            point_up = zeros(numel(vel_up_g)-1,1);
            ss_vel_down = zeros(numel(vel_down_g)-1,4);
            af_vel_down = zeros(numel(vel_down_g)-1,4);
            point_down = zeros(numel(vel_down_g)-1,1);
            
            for i=2:numel(vel_up_g)
                this_ss = ss_all(vel_all>0 & vel_all >=vel_up_g(i-1) & vel_all < vel_up_g(i));
                this_af = af_all(vel_all>0 & vel_all >=vel_up_g(i-1) & vel_all < vel_up_g(i));
                point_up(i-1) = numel(this_ss);
                indx = this_ss<0.5;
                af_vel_up(i-1,1) = nanmean(this_af(indx));
                ss_vel_up(i-1,1) = sum(indx)/numel(this_ss);
                indx = this_ss<1 & this_ss>=0.5;
                af_vel_up(i-1,2) = nanmean(this_af(indx));
                ss_vel_up(i-1,2) = sum(indx)/numel(this_ss);
                indx = this_ss<2 & this_ss>=1;
                af_vel_up(i-1,3) = nanmean(this_af(indx));
                ss_vel_up(i-1,3) = sum(indx)/numel(this_ss);
                indx = this_ss>2;
                af_vel_up(i-1,4) = nanmean(this_af(indx));
                ss_vel_up(i-1,4) = sum(indx)/numel(this_ss); 
            end
            
            for i=2:numel(vel_down_g)
                this_ss = ss_all(vel_all<0 & vel_all <vel_down_g(i-1) & vel_all >= vel_down_g(i));
                this_af = af_all(vel_all<0 & vel_all <vel_down_g(i-1) & vel_all >= vel_down_g(i));
                point_down(i-1) = numel(this_ss);
                indx = this_ss>-0.5;
                af_vel_down(i-1,1) = nanmean(this_af(indx));
                ss_vel_down(i-1,1) = sum(indx)/numel(this_ss);
                indx = this_ss>-1 & this_ss<=-0.5;
                af_vel_down(i-1,2) = nanmean(this_af(indx));
                ss_vel_down(i-1,2) = sum(indx)/numel(this_ss);
                indx = this_ss>-2 & this_ss<=-1;
                af_vel_down(i-1,3) = nanmean(this_af(indx));
                ss_vel_down(i-1,3) = sum(indx)/numel(this_ss);
                indx = this_ss<-2;
                af_vel_down(i-1,4) = nanmean(this_af(indx));
                ss_vel_down(i-1,4) = sum(indx)/numel(this_ss); 
            end
            figure;
            subplot(1,2,1);
            b=barh(ss_vel_up,'stacked');
            labelstr = {'0-2','2-4','4-6','6-8','>8'};
            set(gca,'yticklabel',labelstr); 
            legend('SS <0.5','SS 0.5-1','SS 1-2','SS >2');
            title('Updraft');
            ylabel('W (m/s)');
            
            subplot(1,2,2);
            b=barh(af_vel_up,'group');
            labelstr = {'0-2','2-4','4-6','6-8','>8'};
            set(gca,'yticklabel',labelstr); 
            legend('SS <0.5','SS 0.5-1','SS 1-2','SS >2');
            xlabel('Mean adiabatic fraction');
            ylabel('W (m/s)');
            title('Adia Frac');
            
            figure;
            subplot(1,2,1);
            barh(ss_vel_down,'stacked');
%             labelstr = {'<-8','-8--6','-6--4','-4--2','-2-0'};
            labelstr = {'-2-0','-4--2','-6--4','-8--6','<-8'};
            set(gca,'yticklabel',labelstr)
            legend('SS >-0.5','SS -1--0.5','SS -2--1','SS <-2');
            title('Downdraft');
            ylabel('W (m/s)');
            
            subplot(1,2,2);
            b=barh(af_vel_down,'group');
            labelstr = {'0-2','2-4','4-6','6-8','>8'};
            set(gca,'yticklabel',labelstr); 
            legend('SS <0.5','SS 0.5-1','SS 1-2','SS >2');
            title('Updraft');
            ylabel('Adia Frac');
            xlabel('Mean adiabatic fraction');
       
            
            
            
            % vel vs af vs SS
            af_g = [0.02,0.05,0.1,0.5,1];
            af_up = size(numel(af_g),4);
            af_down = size(numel(af_g),4);
            for i =1:numel(af_g)
                if i==1
                    af_min = 0;
                else
                    af_min = af_g(i-1);
                end
                this_ss = ss_all(vel_all>0 & af_all >=af_min & af_all < af_g(i));
                indx = this_ss<0.5;
                af_up(i,1) = sum(indx)/numel(this_ss);
                indx = this_ss<1 & this_ss>=0.5;
                af_up(i,2) = sum(indx)/numel(this_ss);
                indx = this_ss<2 & this_ss>=1;
                af_up(i,3) = sum(indx)/numel(this_ss);
                indx = this_ss>2;
                af_up(i,4) = sum(indx)/numel(this_ss);
                
                this_ss = ss_all(vel_all<0& af_all >=af_min & af_all < af_g(i));
                indx = this_ss>-0.5;
                af_down(i,1) = sum(indx)/double(numel(this_ss));
                indx = this_ss>-1 & this_ss<=-0.5;
                af_down(i,2) = sum(indx)/double(numel(this_ss));
                indx = this_ss>-2 & this_ss<=-1;
                af_down(i,3) = sum(indx)/double(numel(this_ss));
                indx = this_ss<=-2;
                af_down(i,4) = sum(indx)/double(numel(this_ss));
                 
            end
            
            figure;
            subplot(1,2,1)
            barh(af_up,'stacked');
            set(gca,'yticklabel',{'af <0.02','af 0.02-0.05','af 0.05-0.1','af 0.1-0.5','af 0.5-1','af >1'})
            legend('SS <0.5','SS 0.5-1','SS 1-2','SS >2');
            title('Updraft');
            
            subplot(1,2,2)
            barh(af_down,'stacked');
            set(gca,'yticklabel',{'af <0.02','af 0.02-0.05','af 0.05-0.1','af 0.1-0.5','af 0.5-1','af >1'})
            legend('SS >-0.5','SS -1--0.5','SS -2--1','SS <-2');
            title('Downdraft');
            % vel vs af
            af_g = [0.02,0.05,0.1,0.5,1,10];
            af_up = size(numel(af_g),1);
            af_down = size(numel(af_g),1);
            
            for i =1:numel(af_g)
                if i==1
                    af_min = 0;
                else
                    af_min = af_g(i-1);
                end
                this_af = af_all(vel_all>0);
                indx = this_af >=af_min & this_af < af_g(i);
                af_up(i) = sum(indx)/numel(this_af);
                
                this_af = af_all(vel_all<0);
                indx = this_af >=af_min & this_af < af_g(i);
                af_down(i) = sum(indx)/numel(this_af);
                 
            end
            figure;
            width = 0.4;
            barh([af_up;af_down],'stacked');
            set(gca,'yticklabel',{'updraft','downdraft'})
            legend('<0.02','0.02-0.05','0.05-0.1','0.1-0.5','0.5-1','>1');
            
            % ss vs af
            figure;
            subplot(1,5,1);
            indx = af_all < 0.02;
            nbins = 25;
            histogram(ss_all(indx),nbins);
            xlabel('SS (%)');
            title('AF < 0.0.2');
            
            subplot(1,5,2);
            indx = af_all >= 0.02 & af_all < 0.05;
            histogram(ss_all(indx),nbins);
            xlabel('SS (%)');
            title('AF 0.02-0.05');
            
            subplot(1,5,3);
            indx = af_all >= 0.05 & af_all < 0.1;
            histogram(ss_all(indx),nbins);
            xlabel('SS (%)');
            title('AF 0.05-0.1');
            
            subplot(1,5,4);
            indx = af_all >= 0.1 & af_all < 0.5;
            histogram(ss_all(indx),nbins);
            xlabel('SS (%)');
            title('AF 0.1-0.5');
            
            subplot(1,5,5);
            indx = af_all >= 0.5;
            histogram(ss_all(indx),nbins);
            xlabel('SS (%)');
            title('AF >0.5');
        end
        % same above but to look at the liquid water content
        function explore_lwc_ss()
            opt = 'sim-fixed-dl';
            metfile = '200906160818_1Hz.csv';
            dsdfile = '20090616_30000_1Hz_dsd.csv';
            [SS.d1,nconc.d1,vel.d1,alt.d1,sec.d1,lat.d1,lon.d1,lwc.d1,af.d1,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200906210758_1Hz.csv';
            dsdfile = '20090621_30000_1Hz_dsd.csv';
            [SS.d2,nconc.d2,vel.d2,alt.d2,sec.d2,lat.d2,lon.d2,lwc.d2,af.d2,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200906220753_1Hz.csv';
            dsdfile = '20090622_30000_1Hz_dsd.csv';
            [SS.d3,nconc.d3,vel.d3,alt.d3,sec.d3,lat.d3,lon.d3,lwc.d3,af.d3,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908180730_1Hz.csv';
            dsdfile = '20090818_30000_1Hz_dsd.csv';
            [SS.d4,nconc.d4,vel.d4,alt.d4,sec.d4,lat.d4,lon.d4,lwc.d4,af.d4,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908230850_1Hz.csv';
            dsdfile = '20090823_30000_1Hz_dsd.csv';
            [SS.d5,nconc.d5,vel.d5,alt.d5,sec.d5,lat.d5,lon.d5,lwc.d5,af.d5,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908240858_1Hz.csv';
            dsdfile = '20090824_30000_1Hz_dsd.csv';
            [SS.d6,nconc.d6,vel.d6,alt.d6,sec.d6,lat.d6,lon.d6,lwc.d6,af.d6,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            metfile = '200908250853_1Hz.csv';
            dsdfile = '20090825_30000_1Hz_dsd.csv';
            [SS.d7,nconc.d7,vel.d7,alt.d7,sec.d7,lat.d7,lon.d7,lwc.d7,af.d7,~]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            
            %%%% filter cloud edge points
            day_str = {'d1','d2','d3','d4','d5','d6','d7'};
            for i = 1:numel(day_str)
                clear_sky_indx = lwc.(day_str{i}) < prctile(lwc.(day_str{i}),5);
                clear_sky_sec = sec.(day_str{i})(clear_sky_indx);
                cloud_edge_sec = [];
                for j=1:numel(clear_sky_sec)
                    cloud_edge_sec = horzcat(cloud_edge_sec,clear_sky_sec(j));
                    for k= -2:2
                        if ismember(clear_sky_sec(j)+k, sec.(day_str{i})) && ~ismember(clear_sky_sec(j)+k, cloud_edge_sec)
                            cloud_edge_sec = horzcat(cloud_edge_sec, clear_sky_sec(j)+k);
                        end
                    end
                end
                cloud_edge_indx = zeros(size(sec.(day_str{i})));
                
                for j=1:numel(cloud_edge_sec)
                    this_indx = find(cloud_edge_sec(j)==sec.(day_str{i}));
                    cloud_edge_indx(this_indx) = 1;
                end
                
                ss_cloud_edge.(day_str{i}) = SS.(day_str{i})(logical(cloud_edge_indx));
                ss_cloud_cntr.(day_str{i}) = SS.(day_str{i})(~logical(cloud_edge_indx));
            end
            
            ss_cloud_edge_all = cat(1,ss_cloud_edge.d1,ss_cloud_edge.d2,ss_cloud_edge.d3,ss_cloud_edge.d4,...
                ss_cloud_edge.d5,ss_cloud_edge.d6,ss_cloud_edge.d7);
            ss_cloud_cntr_all = cat(1,ss_cloud_cntr.d1,ss_cloud_cntr.d2,ss_cloud_cntr.d3,ss_cloud_cntr.d4,...
                ss_cloud_cntr.d5,ss_cloud_cntr.d6,ss_cloud_cntr.d7);
            ss_cloud_cntr_all(ss_cloud_cntr_all>20) = nan;
            
            figure;
            hold on;
            histogram(ss_cloud_cntr_all,30,'Normalization','probability');
            histogram(ss_cloud_edge_all,30,'Normalization','probability');
            legend('Non-edge','edge');
            xlabel('SS (%)');
            ylabel('Probability');
            hold off;
            
            ss_all = cat(1,SS.d1,SS.d2,SS.d3,SS.d4,SS.d5,SS.d6,SS.d7);
            af_all = cat(1,af.d1,af.d2,af.d3,af.d4,af.d5,af.d6,af.d7);
            vel_all = cat(1,vel.d1,vel.d2,vel.d3,vel.d4,vel.d5,vel.d6,vel.d7);
            nconc_all = cat(1,nconc.d1,nconc.d2,nconc.d3,nconc.d4,nconc.d5,nconc.d6,nconc.d7);
            alt_all = cat(1,alt.d1,alt.d2,alt.d3,alt.d4,alt.d5,alt.d6,alt.d7);
            lwc_all = cat(1,lwc.d1,lwc.d2,lwc.d3,lwc.d4,lwc.d5,lwc.d6,lwc.d7);
            sec_all = cat(1,sec.d1,sec.d2,sec.d3,sec.d4,sec.d5,sec.d6,sec.d7);
            
            
            %%%% vertical bins
            alt_bins = [2:1:16]*500;
            half_bin = (alt_bins(2)-alt_bins(1))/2;
            lwc_bins = zeros(2,numel(alt_bins));%%% lower than 10th percentile
            ss_bins = zeros(2,numel(alt_bins));
            vel_bins = zeros(2,numel(alt_bins));
            
            for i=1:numel(alt_bins)
                lwc_lower_threshold = prctile(lwc_all(alt_all>alt_bins(i)-half_bin & alt_all<=alt_bins(i)+half_bin),10);
                indx_lower = alt_all>alt_bins(i)-half_bin & alt_all<=alt_bins(i)+half_bin & lwc_all<=lwc_lower_threshold & vel_all>0;
                
                lwc_upper_threshold = prctile(lwc_all(alt_all>alt_bins(i)-half_bin & alt_all<=alt_bins(i)+half_bin),90);
                indx_upper = alt_all>alt_bins(i)-half_bin & alt_all<=alt_bins(i)+half_bin & lwc_all>=lwc_upper_threshold & vel_all>0;
                
                lwc_bins(1,i) = nanmean(lwc_all(indx_lower));
                lwc_bins(2,i) = nanmean(lwc_all(indx_upper));
                ss_bins(1,i) = nanmean(ss_all(indx_lower));
                ss_bins(2,i) = nanmean(ss_all(indx_upper));
                vel_bins(1,i) = nanmean(vel_all(indx_lower));
                vel_bins(2,i) = nanmean(vel_all(indx_upper));
            end
            
            figure;
            subplot(1,3,1);
            hold on;
            line(squeeze(ss_bins(1,:)),alt_bins,'marker','o','color','b');
            line(squeeze(ss_bins(2,:)),alt_bins,'marker','o','color','k');
            xlabel('Supersaturation (%)');
            ylabel('Altitude (m)');
            set(gca,'FontSize', 14);
            hold off;
            
            subplot(1,3,2);
            hold on;
            line(squeeze(lwc_bins(1,:)),alt_bins,'marker','o','color','b');
            line(squeeze(lwc_bins(2,:)),alt_bins,'marker','o','color','k');
            xlabel('Liquid water content (g/m^3)');
            ylabel('Altitude (m)');
            set(gca,'FontSize', 14);
            legend('<10%','>10%')
            hold off;
            
            subplot(1,3,3);
            hold on;
            line(squeeze(vel_bins(1,:)),alt_bins,'marker','o','color','b');
            line(squeeze(vel_bins(2,:)),alt_bins,'marker','o','color','k');
            xlabel('vertical wind  velocity (m/s)');
            ylabel('Altitude (m)');
            set(gca,'FontSize', 14);
            legend('<10%','>10%')
            hold off;
            
        end
        function explore_single_day()
            opt = 'sim-fixed-dl';
            metfile = '200906160818_1Hz.csv';
            dsdfile = '20090616_30000_1Hz_dsd.csv';
            [SS,nconc,vel,alt,sec,lat,lon,lwc,af,rm]=recalculate_CAIPEEX_result.Daily_cal(metfile,dsdfile,opt);
            %%% plot the flight track, determine in/out of cloud
            alt_g = [28:2:70]*100;
            num_g = zeros(size(alt_g));
            fra_g = zeros(size(alt_g));
            for i=1:numel(alt_g)
                indx = alt<alt_g(i)+100 & alt>alt_g(i)-100;
                this_vel = vel(indx);
                num_g(i) = numel(this_vel);
                fra_g(i) = sum(this_vel>0)/num_g(i);
            end
            % figure 3
            ss_ustat = zeros(7,numel(alt_g));
            ss_dstat = zeros(7,numel(alt_g));
            
            for i=1:numel(alt_g)
                indx = alt<alt_g(i)+100 & alt>alt_g(i)-100 & vel>0;
                this_SS = SS(indx);
                ss_ustat(4,i) = numel(this_SS);
                if ss_ustat(4,i)>0
                    ss_ustat(1,i) = min(SS(indx));
                    ss_ustat(2,i) = nanmean(SS(indx));
                    ss_ustat(3,i) = max(SS(indx));
                    ss_ustat(5,i) = sum(this_SS<1)/ss_ustat(4,i);
                    ss_ustat(6,i) = sum(this_SS>=1 & this_SS<2)/ss_ustat(4,i);
                    ss_ustat(7,i) = sum(this_SS>2)/ss_ustat(4,i);
                end
                
                indx = alt<alt_g(i)+100 & alt>alt_g(i)-100 & vel<0;
                this_SS = SS(indx);
                ss_dstat(4,i) = numel(this_SS);
                if ss_dstat(4,i)>0
                    ss_dstat(1,i) = min(SS(indx));
                    ss_dstat(2,i) = nanmean(SS(indx));
                    ss_dstat(3,i) = max(SS(indx));
                    ss_dstat(5,i) = sum(this_SS>-1)/ss_dstat(4,i);
                    ss_dstat(6,i) = sum(this_SS>=-2 & this_SS<=-1)/ss_dstat(4,i);
                    ss_dstat(7,i) = sum(this_SS<-2)/ss_dstat(4,i);
                end
            end
            % figure5
            figure;
            subplot(1,2,1);
            indx = af < 0.1;
            nbins = 25;
            histogram(SS(indx),nbins);
            xlabel('SS (%)');
            title('AF < 0.1');
            
            subplot(1,2,2);
            indx = af>=0.1;
            histogram(SS(indx),nbins);
            xlabel('SS (%)');
            title('AF >= 0.1');
            
            % figure 4
            figure;
            scatter(lwc,af);
            xlabel('LWC (g/m^{3})');
            ylabel('Adia Frac');
            % figure 3
            figure;
            subplot(2,3,1);
            hold on;
            line(squeeze(ss_ustat(1,:)),alt_g,'marker','.','color','b');
            line(squeeze(ss_ustat(2,:)),alt_g,'marker','.','color','k');
            line(squeeze(ss_ustat(3,:)),alt_g,'marker','.','color','r');
            legend('min','mean','max');
            xlabel('SS (%)');
            ylabel('Alt (m)');
            ylim([2.7e3,7.6e3]);
%             title('Updraft');
            hold off;
            
            subplot(2,3,2);
            barh(alt_g,ss_ustat(5:7,:)','stacked');
            legend('<1%','1-2%','>2%');
            ylim([2.7e3,7.6e3]);
            xlabel('Frac');
            
            subplot(2,3,3);
            line(squeeze(ss_ustat(4,:)),alt_g,'marker','o');
            ylim([2.7e3,7.6e3]);
            xlim([0,50]);
            xlabel('# sample points');
            
            subplot(2,3,4);
            hold on;
            line(squeeze(ss_dstat(1,:)),alt_g,'marker','.','color','b');
            line(squeeze(ss_dstat(2,:)),alt_g,'marker','.','color','k');
            line(squeeze(ss_dstat(3,:)),alt_g,'marker','.','color','r');
            legend('min','mean','max');
            xlabel('SS (%)');
            ylabel('Alt (m)');
            ylim([2.7e3,7.6e3]);
%             title('Updraft');
            hold off;
            
            subplot(2,3,5);
            barh(alt_g,ss_dstat(5:7,:)','stacked');
            legend('>-1%','-1--2%','<-2%');
            ylim([2.7e3,7.6e3]);
            xlabel('Frac');
            
            subplot(2,3,6);
            line(squeeze(ss_dstat(4,:)),alt_g,'marker','o');
            ylim([2.7e3,7.6e3]);
            xlabel('# sample points');
            xlim([0,50]);
            
            % figure 2
            indx = alt<6100 & alt>5700;
            sec_s = sec(indx);
            vel_s = vel(indx);
            SS_s = SS(indx);
            nconc_s = nconc(indx);
            lwc_s = lwc(indx);
            af_s = af(indx);
            rm_s = rm(indx);
            
            figure;
            subplot(3,1,1);
            hold on;
            yyaxis left;
            ylabel('W (m/s)');
%             line(1:numel(vel_s),vel_s,'marker','o','color','b','linestyle','-');
            line((sec_s-min(sec_s)),vel_s,'marker','o','color','b','linestyle','none');
            line([1,max(sec_s)],[0,0],'linestyle','--','color','b','linewidth',1);
            yyaxis right;
%             xlim([6.05*60,6.4*60]);
            ylabel('SS (%)');
%             line(1:numel(vel_s),SS_s,'marker','o','color','r','linestyle','-');
            line((sec_s-min(sec_s)),SS_s,'marker','o','color','r','linestyle','none');
            line([1,max(sec_s)],[0,0],'linestyle','--','color','r','linewidth',1);
            hold off;
            
            subplot(3,1,2);
            hold on;
            yyaxis left;
            ylabel('Nc (cm^{-3})');
%             line(1:numel(vel_s),nconc_s,'marker','o','color','b','linestyle','-');
            line((sec_s-min(sec_s)),nconc_s,'marker','o','color','b','linestyle','none');
            line([1,max(sec_s)],[mean(nconc_s),mean(nconc_s)],'linestyle','--','color','b','linewidth',1);
            yyaxis right;
            ylabel('Mean radius (microns)');
            ylim([0,20]);
%             xlim([6.05*60,6.4*60]);
%             line(1:numel(vel_s),rm_s,'marker','o','color','r','linestyle','-');
            line((sec_s-min(sec_s)),rm_s,'marker','o','color','r','linestyle','none');
            line([1,max(sec_s)],[mean(rm_s),mean(rm_s)],'linestyle','--','color','r','linewidth',1);
            hold off;
            
            subplot(3,1,3);
            xlabel('sample points')
            hold on;
            yyaxis left;
            ylabel('LWC (g/m^{3})');
            ylim([0,4]);
%             line(1:numel(vel_s),lwc_s,'marker','o','color','b','linestyle','-');
            line((sec_s-min(sec_s)),lwc_s,'marker','o','color','b','linestyle','none');
            line([1,max(sec_s)],[0,0],'linestyle','--','color','b','linewidth',1);
            yyaxis right;
            ylabel('Adiat Frac');
%             xlim([6.05*60,6.4*60]);
            
%             line(1:numel(vel_s),af_s,'marker','o','color','r','linestyle','-');
            line((sec_s-min(sec_s)),af_s,'marker','o','color','r','linestyle','none');
            line([1,max(sec_s)],[0.02,0.02],'linestyle','--','color','r','linewidth',1);
%             xlim([6.05*60,6.4*60]);
            hold off;
            
            
            % figure 1
            figure;
            subplot(1,3,1);
            scatter(vel,alt,[],(sec-min(sec))/60,'filled');
            h=colorbar;
            ylabel(h, 'Sample time (min)')
            xlabel('W (m/s)');
            ylabel('Altitude (m)');
            line([0,0],[0,7e3],'linestyle','--');
            ylim([2.7e3,7.6e3]);
            
            subplot(1,3,2);
            line(num_g,alt_g,'marker','o');
            ylim([2.7e3,7.6e3]);
            xlabel('# sample points');
%             ylabel('Altitude (m)');
            
            subplot(1,3,3);
            barh(alt_g,[fra_g;1-fra_g]','stacked');
            ylim([2.7e3,7.6e3]);
            legend('Updraft','Downdraft');
            line([0.5,0.5],[2e3,8e3],'color','k','linestyle','--','linewidth',3)
            
            figure;
            hold on;
            scatter3(lon,lat,alt,[],sec/3600,'filled');
            h = surface([lon(:), lon(:)], [lat(:), lat(:)], [alt(:), alt(:)], ...
                [sec(:)/3600, sec(:)/3600], 'EdgeColor','flat', 'FaceColor','none');
            
            
%             c = 1:numel(t);      %# colors
%             h = surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
%                 [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none');
%             colormap( jet(numel(t)) )
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
    end
    
    
end

