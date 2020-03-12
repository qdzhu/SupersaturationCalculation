classdef explore_model_results
    properties(Constant = true)
    end
    methods(Static)
        
        function read_wrf_var()
            filepath = '/Volumes/share2/USERS/ZhuQ/Cloud microphysics/C_BG/';
            filedir = dir(fullfile(filepath,'wrfout_subset_d01_2014-03-17_14*'));
            wrfinfo = ncinfo(fullfile(filepath,filedir(1).name));
            lon = ncread(wrfinfo.Filename,'XLONG');
            lat = ncread(wrfinfo.Filename,'XLAT');
            ss = ncread(wrfinfo.Filename,'SSW');
            vel = ncread(wrfinfo.Filename,'W');
            nc = ncread(wrfinfo.Filename,'QNCLOUD');
            tke = ncread(wrfinfo.Filename,'TKE');
            z = read_wrf_preproc(wrfinfo.Filename,'z');
            save('wrf_pi_14_bg','lon','lat','ss','vel','nc','tke','z');
            
        end
        
        function test_steady_state_ss()
            
            [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=recalculate_CAIPEEX_result.Constant; 
            
            filepath = '/Volumes/share2/USERS/ZhuQ/Cloud microphysics/';
            wrfgeos = ncinfo(fullfile(filepath,'wrf_geosinfo'));
            lon = ncread(wrfgeos.Filename,'XLONG');
            lon = squeeze(lon(:,:,1));
            lat = ncread(wrfgeos.Filename,'XLAT');
            lat = squeeze(lat(:,:,1));
            height = read_wrf_preproc(wrfgeos.Filename,'z_center');
            height = squeeze(height(:,:,:,1));
            wrf_T = ncread(wrfgeos.Filename,'T');
            
            wrfinfo = ncinfo(fullfile(filepath, 'C_PI','wrfout_subset_d01_2014-03-17_18'));
            wrf_P = ncread(wrfinfo.Filename, 'P');
            wrf_PB = ncread(wrfinfo.Filename,'PB');
            
            temp = convert_wrf_temperature(wrf_T, wrf_P, wrf_PB);
            temp = squeeze(temp(:,:,:,6)); %%% we only select one time step for testing purpose
            pres = read_wrf_preproc(wrfinfo.Filename,'pres'); %%% hpa
            pres = squeeze(pres(:,:,:,6));
            
            qncloud = ncread(wrfinfo.Filename, 'QNCLOUD');
            qncloud = squeeze(qncloud(:,:,:,6));
            
            bin_diameter = zeros(17,1);
            bin_diameter(1) = 2;
            for i=2:15
                bin_diameter(i) = bin_diameter(i-1)*2^(1/3);
            end
            
            varname = num2str([01:33].','ff1i%02d');
            total_weighted_var = zeros(size(qncloud));
            total_var = zeros(size(qncloud));
            for i=1:15
                this_sd = ncread(wrfinfo.Filename,varname(i,:));
                total_weighted_var = total_weighted_var + squeeze(this_sd(:,:,:,6))/2^(i-1)*bin_diameter(i);
                total_var = total_var + squeeze(this_sd(:,:,:,6))/2^(i-1);
            end
            aver_rccn = total_weighted_var./total_var;
            save('test_steady_state_ss','temp','pres','qncloud','aver_rccn','bin_diameter');
        end
        
        function test_steady_state_ss_p2()
            filepath = '/Volumes/share2/USERS/ZhuQ/Cloud microphysics/';
            wrfgeos = ncinfo(fullfile(filepath,'wrf_geosinfo'));
            lon = ncread(wrfgeos.Filename,'XLONG');
            lon = squeeze(lon(:,:,1));
            lat = ncread(wrfgeos.Filename,'XLAT');
            lat = squeeze(lat(:,:,1));
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            %%%% find the index to filter the study domain
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            height = read_wrf_preproc(wrfgeos.Filename,'z_center');
            height = squeeze(height(:,:,:,1));
            height = permute(height,[3,1,2]);
            height = height(:,regionindx);
            
            wrfinfo = ncinfo(fullfile(filepath, 'C_PI','wrfout_subset_d01_2014-03-17_18'));
            vel = ncread(wrfinfo.Filename, 'W');
            vel = squeeze(vel(:,:,:,6));
            vel = (vel(:,:,1:end-1)+vel(:,:,2:end))/2;
            vel = permute(vel,[3,1,2]);
            vel = vel(:,regionindx);
            SSW = ncread(wrfinfo.Filename,'SSW');
            SSW = squeeze(SSW(:,:,:,6));
            SSW = permute(SSW,[3,1,2]);
            SSW = SSW(:,regionindx);
            
            [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=recalculate_CAIPEEX_result.Constant; 
            
            data = load('test_steady_state_ss');
            temp = data.temp;%%k
            temp = permute(temp,[3,1,2]);
            temp = temp(:,regionindx);
            
            pres = data.pres;%hpa
            pres = permute(pres,[3,1,2]);
            pres = pres(:,regionindx);
            
            qncloud = data.qncloud; %%%% kg-1
            qncloud = permute(qncloud,[3,1,2]);
            qncloud = qncloud(:,regionindx);
            
            Tc = temp-To; % celsius degree
            aver_rccn = data.aver_rccn;%%%% micron
            aver_rccn = permute(aver_rccn,[3,1,2]);
            aver_rccn = aver_rccn(:,regionindx);
            
            ndens = pres*100./(Ra*temp);
            qncloud = qncloud.*ndens/1e6;
            
            L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
            D=(2.26e-5+1.5e-7*Tc)*Po./(pres*1e2);  % diffusion coeff.
            A = (g*L./(Cpa*Rv*temp.^2)-g./(Ra*temp));
            SS = A.*vel./(4*pi*D.*aver_rccn.*qncloud);% quasi-state supersaturation
            
            figure;
            indx = qncloud > prctile(qncloud(:),99);%% 2.92
            scatter(SSW(indx)*100, SS(indx)*100,10,vel(indx));
            h = colorbar;
            ylabel(h, 'W (m/s)');
            plot_fit_line(SSW(indx)*100, SS(indx)*100);
            xlabel('SS from WRF (%) ');
            ylabel('Steady-state SS (%)');
            ylim([-10,15]);
            set(gca, 'fontsize', 14)
            
            figure;
            scatter(SSW(indx)*100, SS(indx)*100,10,qncloud(indx));
            h = colorbar;
            ylabel(h, 'Droplet number concentration(cm^{-3})');
            plot_fit_line(SSW(indx)*100, SS(indx)*100);
            xlabel('SS from WRF (%) ');
            ylabel('Steady-state SS (%)');
            ylim([-10,15]);
            set(gca, 'fontsize', 14)
            
        end
        
        function read_wrf_study_domain()
            filepath = '/Volumes/share2/USERS/ZhuQ/Cloud microphysics/';
            %%%%% read the longitude, latitude, height from wrf
            wrfgeos = ncinfo(fullfile(filepath,'wrf_geosinfo'));
            lon = ncread(wrfgeos.Filename,'XLONG');
            lon = squeeze(lon(:,:,1));
            lat = ncread(wrfgeos.Filename,'XLAT');
            lat = squeeze(lat(:,:,1));
            height = read_wrf_preproc(wrfgeos.Filename,'z_center');
            height = squeeze(height(:,:,:,1));
            
            modelopt = {'C_PI','C_BG'};
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            %%%% find the index to filter the study domain
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            
            study_hour = 14:19;

            
%             data_fields = {'q','temp'};
%             data = make_empty_struct_from_cell(data_fields,nan(numel(study_hour)*12,sum(regionindx(:)) ,size(height,3)));
            varname = num2str([01:33].','ff1i%02d');
            
            bin_diameter = zeros(15,1);
            bin_diameter(1) = 2;
            for i=2:15
                bin_diameter(i) = bin_diameter(i-1)*2^(1/3);
            end
            
            for j=1:2
                for i=1:numel(study_hour)
                    filedir = dir(fullfile(filepath,modelopt{j},sprintf('wrfout_subset_d01_2014-03-17_%02d*', study_hour(i))));
                    wrfinfo = ncinfo(fullfile(filepath,modelopt{j},filedir(1).name));
                    for k=1:15
                        display(k);
                        thisvarname = varname(k,:);
                        this_sd = ncread(wrfinfo.Filename,varname(k,:))/2^(k-1);
                        data.(thisvarname)(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(this_sd,regionindx);
                    end
                end
                if j ==1
                    PI = data;
                else
                    BG = data;
                end
                
            end
            save('read_wrf_study_domain_dsd','BG','PI','lon','lat','height','bin_diameter');
            
            for j=1:2
                for i=1:numel(study_hour)
                    filedir = dir(fullfile(filepath,modelopt{j},sprintf('wrfout_subset_d01_2014-03-17_%02d*', study_hour(i))));
                    wrfinfo = ncinfo(fullfile(filepath,modelopt{j},filedir(1).name));
                    q = ncread(wrfinfo.Filename,'QVAPOR');
                    data.q(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(q, regionindx);
                    
                    filedir = dir(fullfile(filepath,modelopt{j},sprintf('wrfout_subset_T_d01_2014-03-17_%02d*', study_hour(i))));
                    wrfinfo = ncinfo(fullfile(filepath,modelopt{j},filedir(1).name));
                    temp = ncread(wrfinfo.Filename,'T');
                    data.temp(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(temp,regionindx);
                end
                 if j ==1
                    PI = data;
                else
                    BG = data;
                end
            end
            save('read_extra_wrf_study_domain','BG','PI','lon','lat','height');
            
            for j =1:2
                for i=1:numel(study_hour)
                    filedir = dir(fullfile(filepath,modelopt{j},sprintf('wrfout_subset_d01_2014-03-17_%02d*', study_hour(i))));
                    wrfinfo = ncinfo(fullfile(filepath,modelopt{j},filedir(1).name));
                    ss = ncread(wrfinfo.Filename,'SSW');
                    data.ss(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(ss, regionindx);
                    ss_nconc = ncread(wrfinfo.Filename,'SSW_RATIO');
                    data.ss_nconc(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(ss_nconc, regionindx);
                    nconc = ncread(wrfinfo.Filename, 'QNCLOUD');
                    data.nconc(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(nconc,regionindx);
                    w = ncread(wrfinfo.Filename, 'W');
                    w_cut = explore_model_results.cut_by_region_indx(w,regionindx);
                    data.w(12*i-11:12*i,:,:) = (w_cut(:,:,1:end-1)+w_cut(:,:,2:end))/2;
                    latentrate = ncread(wrfinfo.Filename,'TEMPDIFFL');
                    data.lc_water(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(latentrate,regionindx);
                    filedir = dir(fullfile(filepath,modelopt{j},sprintf('wrfout_add_subset_d01_2014-03-17_%02d*', study_hour(i))));
                    wrfinfo = ncinfo(fullfile(filepath,modelopt{j},filedir(1).name));
                    latentrate_ice = ncread(wrfinfo.Filename,'TEMPDIFFI');
                    data.lc_ice(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(latentrate_ice,regionindx);
                    latentrate_freeze = ncread(wrfinfo.Filename,'TFRZMELT');
                    data.lc_freeze(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(latentrate_freeze,regionindx);
                    latentrate_rim = ncread(wrfinfo.Filename,'TEMPRIM');
                    data.lc_rim(12*i-11:12*i,:,:) = explore_model_results.cut_by_region_indx(latentrate_rim,regionindx);
                    
                end  
                if j ==1
                    PI = data;
                else
                    BG = data;
                end
            end
            save('read_wrf_study_domain_dsd','BG','PI','lon','lat','height');
        end
        
        function make_aver_radius()
            % data struct to save the output
            output = load('wrf_supersaturation_comparison');
            % read in the did data
            data = load('read_wrf_study_domain_dsd');
            fieldnames = {'PI','BG'};
            varname = num2str([01:33].','ff1i%02d');
            bin_diameter = data.bin_diameter;
            for i = 1:numel(fieldnames)
                strc = data.(fieldnames{i});
                total_weight = zeros(size(strc.(varname(1,:))));
                total_var = zeros(size(total_weight));
                for k = 1:15
                    total_weight = total_weight + strc.(varname(k,:))*bin_diameter(k);
                    total_var = total_var + strc.(varname(k,:));
                end
                aver_radius = 0.5* total_weight./total_var;
                output.(fieldnames{i}).aver_radius = aver_radius;
                output.(fieldnames{i}).S_qss = output.(fieldnames{i}).S_qss_times_rmean./aver_radius;
            end
            
        end
        
        function make_qss_ss()
            output = load('wrf_supersaturation_comparison');
            output = output.output;
            data = load('read_wrf_study_domain_fan');
            extra = load('read_wrf_study_domain');
            pres = data.pres;
            pres = repmat(pres,6,1);
            [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=recalculate_CAIPEEX_result.Constant; 
            
            fieldnames = {'PI','BG'};
            for i=1:numel(fieldnames)
                temp = data.(fieldnames{i}).temp;
                ndens = pres*100./(Ra*temp);
                Tc = temp - To;
                L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
                D=(2.26e-5+1.5e-7*Tc)*Po./(pres*1e2);  % diffusion coeff.
                A = (g*L./(Cpa*Rv*temp.^2)-g./(Ra*temp));
                vel = extra.(fieldnames{i}).w;
                aver_diameter = output.(fieldnames{i}).aver_radius*2;
                qncloud = extra.(fieldnames{i}).nconc.*ndens/1e6;
                SS = A.*vel./(4*pi*D.*aver_diameter.*qncloud);% quasi-state supersaturation
                output.(fieldnames{i}).S_qss = SS;
                output.(fieldnames{i}).qncloud = qncloud;
                output.(fieldnames{i}).vel = vel;
            end
            
            save('wrf_supersaturation_comparison_v1.mat','output');
        end
        
        function make_fan_ss()
            data = load('read_extra_wrf_study_domain.mat');
            height = data.height;
            lon = data.lon;
            lat = data.lat;
            altitude = zeros(size(height,3),1);
            PI = data.PI;
            BG = data.BG;
            %%%% find the index to filter the study domain
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            
            filepath = '/Volumes/share2/USERS/ZhuQ/Cloud microphysics/';
            %%%%% read the longitude, latitude, height from wrf
%             wrfgeos = ncinfo(fullfile(filepath,'wrf_geosinfo'));
%             lon = ncread(wrfgeos.Filename,'XLONG');
%             lon = squeeze(lon(:,:,1));
%             lat = ncread(wrfgeos.Filename,'XLAT');
%             lat = squeeze(lat(:,:,1));
            
%             height = read_wrf_preproc(wrfgeos.Filename,'z_center');
%             height = squeeze(height(:,:,:,1));
            
            pres = read_wrf_preproc(fullfile(filepath,'wrfout_presinfo'),'pres');
            pres = explore_model_results.cut_by_region_indx(pres,regionindx);
            P = ncread(fullfile(filepath,'wrfout_presinfo'),'P');
            PB = ncread(fullfile(filepath,'wrfout_presinfo'),'PB');
            P = explore_model_results.cut_by_region_indx(P,regionindx);
            PB = explore_model_results.cut_by_region_indx(PB,regionindx);
            
            PI.ss_fan = nan(size(PI.temp));
            BG.ss_fan = nan(size(BG.temp));
            for i=1:6
                fprintf('Implements %d\n',i);
                this_temp = convert_wrf_temperature(PI.temp(12*i-11:12*i,:,:), P,PB);
                this_es=6.112*exp(17.67*(this_temp-273.15)./(this_temp-29.65));
                this_qs = 0.622*this_es./(pres-this_es);
                this_q = PI.q(12*i-11:12*i,:,:);
                this_ss = this_q./this_qs-1;
                PI.ss_fan(12*i-11:12*i,:,:) = this_ss;
                PI.temp(12*i-11:12*i,:,:) = this_temp;
                
                this_temp = convert_wrf_temperature(BG.temp(12*i-11:12*i,:,:), P,PB);
                this_es=6.112*exp(17.67*(this_temp-273.15)./(this_temp-29.65));
                this_qs = 0.622*this_es./(pres-this_es);
                this_q = BG.q(12*i-11:12*i,:,:);
                this_ss = this_q./this_qs-1;
                BG.ss_fan(12*i-11:12*i,:,:) = this_ss;
                BG.temp(12*i-11:12*i,:,:) = this_temp;
            end
            
            save('read_wrf_study_domain_fan','BG','PI','lon','lat','height','pres');
        end
        % when we compare the fan's ss against wrf ss, we find out fan's ss
        % scatter wildly when wrf gives zeros values. The method here is to
        % figure out why
        function match_fan_wrf_ss()
            data = load('read_wrf_study_domain.mat');
            wrf = data.PI;
            wrf_ss = data.PI.ss;
            wrf_ss_nconc = data.PI.ss_nconc;
            
            data = load('read_wrf_study_domain_fan.mat');
            pres = data.pres;
            pres = repmat(pres,6,1,1);
            PI = data.PI;
            fan_ss = PI.ss_fan;
            
            % wrf_ss against wrf_ss_nconc
            indx = wrf_ss_nconc ==0 & wrf_ss ~= 0;
            wpres = pres(indx);
            wtemp = PI.temp(indx);
            
            indx = wrf_ss == 0 | wrf_ss_nconc ==0;
            wrf_ss(indx) = nan;
            wrf_ss_nconc(indx) = nan;
            % first look at the datasets
            figure;
            subplot(2,2,1);
            hold on;
            scatter(wrf_ss(:), wrf_ss_nconc(:),15, 'b', 'filled');
            scatter(wrf_ss(indx), wrf_ss_nconc(indx),10,'y');
            line([-1,0.8],[-1,0.8],'color','r','linewidth',1);
            xlabel('WRF SS_1');
            ylabel('WRF SS_2');
            hold off;
            
            subplot(2,2,2);
            scatter(fan_ss(:), wrf_ss(:),15, 'b', 'filled');
            line([-1,0.8],[-1,0.8],'color','r','linewidth',1);
            ylabel('WRF SS_1');
            xlabel('Fan SS');
            
            subplot(2,2,3);
            scatter(fan_ss(:), wrf_ss_nconc(:),15, 'b', 'filled');
            line([-1,0.8],[-1,0.8],'color','r','linewidth',1);
            ylabel('WRF SS_2');
            xlabel('Fan SS');
            
        end
        
        function value_2d = reconst_spatial_pattern(value, lon, lat)
            
            modelopt = {'C_PI','C_BG'};
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            %%%% find the index to filter the study domain
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            s_lon = lon(regionindx);
            s_lat = lat(regionindx);
             
            value_2d = nan(size(lon));
            value_2d(regionindx) = value;
        end
        
        function make_cmp_3d()
            % get longitude and latitude info
            data = load('read_wrf_study_domain_fan.mat');
            lon = data.lon;
            lat = data.lat;
            % extract out the datasets
%             data = load('wrf_supersaturation_comparison_v1.mat');
            output = load('wrf_lwc_data.mat');
            bg_lwc = squeeze(output.wrf_lwc_data.BG_lwc);
            pi_lwc = squeeze(output.wrf_lwc_data.PI_lwc);
%             data.output.PI.lwc = pi_lwc;
%             data.output.BG.lwc = bg_lwc;
            
            fieldnames = {'temp'};
            PI = make_empty_struct_from_cell(fieldnames);
            for i=1:numel(fieldnames)
                PI.(fieldnames{i}) = nan(size(bg_lwc,1),size(lon,1),size(lon,2),size(bg_lwc,3));
                for j=1:size(bg_lwc,1)
                    for k=1:size(bg_lwc,3)
                        PI.(fieldnames{i})(j,:,:,k) = ...
                            explore_model_results.reconst_spatial_pattern...
                            (data.PI.(fieldnames{i})(j,:,k), lon, lat);
                    end
                end
            end
            save('read_wrf_study_domain_3d_temp_pi.mat','PI','-v7.3');
            for i=1:numel(fieldnames)
                PI.(fieldnames{i}) = nan(size(bg_lwc,1),size(lon,1),size(lon,2),size(bg_lwc,3));
                for j=1:size(bg_lwc,1)
                    for k=1:size(bg_lwc,3)
                        PI.(fieldnames{i})(j,:,:,k) = ...
                            explore_model_results.reconst_spatial_pattern...
                            (data.output.PI.(fieldnames{i})(j,:,k), lon, lat);
                    end
                end
            end
            save('read_wrf_study_domain_3d_pi.mat','PI','-v7.3');
            for i=1:numel(fieldnames)
                BG.(fieldnames{i}) = nan(size(bg_lwc,1),size(lon,1),size(lon,2),size(bg_lwc,3));
                for j=1:size(bg_lwc,1)
                    for k=1:size(bg_lwc,3)
                        BG.(fieldnames{i})(j,:,:,k) = ...
                            explore_model_results.reconst_spatial_pattern...
                            (data.output.BG.(fieldnames{i})(j,:,k), lon, lat);
                    end
                end
            end
            % lwc
            
            save('read_wrf_study_domain_3d_bg.mat','BG','-v7.3');
            
        end
        
        % plot procedure
        function plot_cmp_3d()
            data = load('read_wrf_study_domain_fan.mat');
            lon = data.lon;
            lat = data.lat;
            height = data.height;
            
            data = load('read_wrf_study_domain_3d_pi');
            pi = data.PI;
         
            vel = pi.vel;
            vel(pi.lwc <= 1e-5) = nan;
            s_fan = pi.S_fan;
            s_fan(pi.lwc <= 1e-5) = nan;
            s_qss = pi.S_qss;
            s_qss(pi.lwc <= 1e-5) = nan;
            
            data = load('read_wrf_study_domain_3d_temp_pi');
            temp = data.PI.temp;
            temp(pi.lwc <= 1e-5) = nan;
            
            lon_bdy = [-60.9, -60.4];
            lat_bdy = [-3.4, -3.0];
            
            % relationship between temperature and s_fan
            
            writerObj = VideoWriter('PI_filter_xz.avi');
            writerObj.FrameRate = 4;
            open(writerObj);

            for i=1:size(vel,1)
                
                xlon = repmat(lon,1,1,66);
                xlat = repmat(lat,1,1,66);
                figure(1);
                
                subplot(1,3,1);
                xvel = squeeze(vel(i,:,:,:));
                scatter3(xlon(:), xlat(:), height(:)/1e3,  5, xvel(:));
                h = colorbar;
                caxis([-5,5]);
                colormap(blue_red_cmap);
                xlim(lon_bdy);
                ylim(lat_bdy);
                zlim([0,11]);
                xlabel('Longitude');
                ylabel('Latitude');
                zlabel('Height (km)');
                ylabel(h,'W (m/s)');
                title('Vertical wind velocity (m/s)');
                view(90,0);
                
                subplot(1,3,2);
                xs_fan = squeeze(s_fan(i,:,:,:));
                scatter3(xlon(:), xlat(:), height(:)/1e3,  5, xs_fan(:)*100);
                h = colorbar;
                caxis([-20,20]);
                colormap(blue_red_cmap);
                xlim(lon_bdy);
                ylim(lat_bdy);
                zlim([0,11]);
                xlabel('Longitude');
                ylabel('Latitude');
                zlabel('Height (km)');
                ylabel(h,'Supersaturation (%)');
                title('Supersaturation from fan calculation (%)');
                view(90,0);
                
                subplot(1,3,3);
                xs_qss = squeeze(s_qss(i,:,:,:));
                scatter3(xlon(:), xlat(:), height(:)/1e3,  5, xs_qss(:)*100);
                h = colorbar;
                caxis([-20,20]);
                colormap(blue_red_cmap);
                xlim(lon_bdy);
                ylim(lat_bdy);
                zlim([0,11]);
                xlabel('Longitude');
                ylabel('Latitude');
                zlabel('Height (km)');
                ylabel(h,'Supersaturation (%)');
                title('Supersaturation from qss calculation (%)');
                view(90,0);
                
                set(gcf,'Position',[100 100 1500 400]);
                frame = getframe(gcf);
                writeVideo(writerObj, frame);
            end
            
            close(writerObj);
        end
        
        function plot_cmp_temp_3d()
            data = load('read_wrf_study_domain_fan.mat');
            lon = data.lon;
            lat = data.lat;
            height = data.height;
            
            data = load('read_wrf_study_domain_3d_pi');
            pi = data.PI;
         
            vel = pi.vel;
            vel(pi.lwc <= 1e-5) = nan;
            s_fan = pi.S_fan;
            s_fan(pi.lwc <= 1e-5) = nan;
            s_qss = pi.S_qss;
            s_qss(pi.lwc <= 1e-5) = nan;
            
            data = load('read_wrf_study_domain_3d_temp_pi');
            temp = data.PI.temp;
            temp(pi.lwc <= 1e-5) = nan;
            
            lon_bdy = [-60.9, -60.4];
            lat_bdy = [-3.4, -3.0];
            
            % relationship between temperature and s_fan
            
            writerObj = VideoWriter('PI_temp_ssfan.avi');
            writerObj.FrameRate = 4;
            open(writerObj);

            for i=1:size(vel,1)
                
                xlon = repmat(lon,1,1,66);
                xlat = repmat(lat,1,1,66);
                figure(1);
                
                subplot(1,3,1);
                xtemp = squeeze(temp(i,:,:,:));
                scatter3(xlon(:), xlat(:), height(:)/1e3,  5, xtemp(:)-273);
                h = colorbar;
                caxis([-20,20]);
                colormap(blue_red_cmap);
                xlim(lon_bdy);
                ylim(lat_bdy);
                zlim([0,11]);
                xlabel('Longitude');
                ylabel('Latitude');
                zlabel('Height (km)');
                ylabel(h,'Temperature - 273 (K)');
                title('Temperature -273 (K)');
          
                
                subplot(1,3,2);
                xs_fan = squeeze(s_fan(i,:,:,:));
                scatter3(xlon(:), xlat(:), height(:)/1e3,  5, xs_fan(:)*100);
                h = colorbar;
                caxis([-20,20]);
                colormap(blue_red_cmap);
                xlim(lon_bdy);
                ylim(lat_bdy);
                zlim([0,11]);
                xlabel('Longitude');
                ylabel('Latitude');
                zlabel('Height (km)');
                ylabel(h,'Supersaturation (%)');
                title('Supersaturation from fan calculation (%)');
            
                
                subplot(1,3,3);
                xs_qss = squeeze(s_qss(i,:,:,:));
                scatter3(xlon(:), xlat(:), height(:)/1e3,  5, xs_qss(:)*100);
                h = colorbar;
                caxis([-20,20]);
                colormap(blue_red_cmap);
                xlim(lon_bdy);
                ylim(lat_bdy);
                zlim([0,11]);
                xlabel('Longitude');
                ylabel('Latitude');
                zlabel('Height (km)');
                ylabel(h,'Supersaturation (%)');
                title('Supersaturation from qss calculation (%)');
                
                set(gcf,'Position',[100 100 1500 400]);
                frame = getframe(gcf);
                writeVideo(writerObj, frame);
            end
            
            close(writerObj);
        end
        
        function plot_cmp_fan_qss()
            data = load('wrf_supersaturation_comparison_v1.mat');
            PI = data.output.PI;
            
            indx = PI.qncloud(:) > 10; % cloud droplet number concentration>10 particles/cm^3
            figure;
            subplot(1,2,1);
            hold on;
            scatter(PI.S_fan(indx)*100, PI.S_qss(indx)*100,3,'MarkerFaceColor','none',...
                'LineWidth',0.2,'MarkerEdgeColor','blue');
            line([-20,20],[-20,20],'linestyle','-','color','r','linewidth',2);
            line([-20,20],[0,0],'linestyle','--','color','r','linewidth',1);
            line([0,0],[-20,80],'linestyle','--','color','r','linewidth',1);
            xlim([-15,15]);
            xlabel('Supersaturatin from Fan (%)');
            ylabel('Supersaturation from Steady State (%)');
            hold off;
            
            subplot(1,2,2);
            hold on;
            indx = PI.qncloud(:) > 3 & PI.S_qss(:) >0 & PI.S_fan(:)>0;
            scatter(PI.S_fan(indx)*100, PI.S_qss(indx)*100,3,'MarkerFaceColor','none',...
                'LineWidth',0.2,'MarkerEdgeColor','blue');
            plot_fit_line(PI.S_fan(indx)*100, PI.S_qss(indx)*100);
            xlabel('Supersaturatin from Fan (%)');
            ylabel('Supersaturation from Steady State (%)');
            hold off;
        end
        
        function plot_cmp_fan_qss_v1()
            data = load('wrf_supersaturation_comparison_v1.mat');
            PI = data.output.PI;
            
            output = load('wrf_lwc_data.mat');
            lwc = squeeze(output.wrf_lwc_data.PI_lwc);
            indx = lwc(:) > 1e-5 & PI.qncloud(:) > 3 & PI.vel(:)>2; % cloud droplet number concentration>10 particles/cm^3
            figure;
            subplot(1,2,1);
            hold on;
            scatter(PI.S_fan(indx)*100, PI.S_qss(indx)*100,3,'MarkerFaceColor','none',...
                'LineWidth',0.2,'MarkerEdgeColor','blue');
            line([-20,20],[-20,20],'linestyle','-','color','r','linewidth',2);
            line([-20,20],[0,0],'linestyle','--','color','r','linewidth',1);
            line([0,0],[-20,80],'linestyle','--','color','r','linewidth',1);
            ylim([-15,15]);
            xlabel('Supersaturatin from Fan (%)');
            ylabel('Supersaturation from Steady State (%)');
            hold off;
            
            subplot(1,2,2);
            hold on;
            scatter(PI.S_fan(indx)*100, PI.S_qss(indx)*100,3,lwc(indx));
            line([-20,20],[-20,20],'linestyle','-','color','r','linewidth',2);
            line([-20,20],[0,0],'linestyle','--','color','r','linewidth',1);
            line([0,0],[-20,80],'linestyle','--','color','r','linewidth',1);
            ylim([-15,15]);
            xlabel('Supersaturatin from Fan (%)');
            ylabel('Supersaturation from Steady State (%)');
            hold off;
            
            subplot(1,2,2);
            hold on;
            indx = PI.qncloud(:) > 3 & PI.S_qss(:) >0 & PI.S_fan(:)>0;
            scatter(PI.S_fan(indx)*100, PI.S_qss(indx)*100,3,'MarkerFaceColor','none',...
                'LineWidth',0.2,'MarkerEdgeColor','blue');
            plot_fit_line(PI.S_fan(indx)*100, PI.S_qss(indx)*100);
            xlabel('Supersaturatin from Fan (%)');
            ylabel('Supersaturation from Steady State (%)');
            hold off;
        end
        
        function plot_wrf_study_domain()
            data = load('read_wrf_study_domain_fan.mat');
            height = data.height;
            lon = data.lon;
            lat = data.lat;
            altitude = zeros(size(height,3),1);
            PI = data.PI;
            BG = data.BG;
            %%%% find the index to filter the study domain
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            
            for i=1:size(height,3)
                height_cut = squeeze(height(:,:,i));
                altitude(i) = nanmean(height_cut(regionindx));
            end
            %%%% top 10 percentiles for the updrafts with w>2 m/s during
            %%%% 1400-1900 UTC from the convective clouds
            opt = 'alltime';
            uv = explore_model_results.filter_fan(PI, BG, 'w', altitude,opt);
            ss = explore_model_results.filter_fan(PI,BG, 'ss', altitude,opt);
            ss_nconc = explore_model_results.filter_fan(PI,BG,'ss_nconc', altitude,opt);
            ss_fan = explore_model_results.filter_fan(PI,BG,'ss_fan',altitude,opt);
%             ss_weighted = explore_model_results.filter_fan(PI, BG, 'weighted_ss', altitude,opt);
            %compare ss from pan calculation and from model output
            
            figure;
            hold on;
            scatter(PI.ss(:)*100,PI.ss_fan(:)*100);
            xlabel('SS (%) from WRF');
            ylabel('SS (%) from Fan calculation');
            line([-50,50],[-50,50],'linestyle','--','color','r');
            hold off;
            
            figure;
            subplot(2,2,1);
            hold on;
            line(uv.pi, uv.alt/1e3,'Linestyle','--', 'color','b','linewidth',2);
            line(uv.bg, uv.alt/1e3,'Linestyle','-', 'color','b','linewidth',2);
            xlabel('Updraft velocity (m/s)');
            ylabel('Height (km)');
            legend('C_ PI','C_ BG');
            ylim([0,15]);
            xlim([0,10]);
            hold off;
            
            subplot(2,2,2);
            hold on;
            line(ss.pi*100, ss.alt/1e3,'Linestyle','--', 'color','b','linewidth',2);
            line(ss.bg*100, ss.alt/1e3,'Linestyle','-', 'color','b','linewidth',2);
            xlabel('Supersaturation (%)');
            ylabel('Height (km)');
            legend('C_ PI','C_ BG');
%             xlim([0,20]);
            ylim([0,15]);
            xlim([0,20]);
            title('SS from WRF output');
            hold off;
            
            subplot(2,2,3);
            hold on;
            line(ss_nconc.pi*100, ss_fan.alt/1e3,'Linestyle','--', 'color','b','linewidth',2);
            line(ss_nconc.bg*100, ss_nconc.alt/1e3,'Linestyle','-', 'color','b','linewidth',2);
            xlabel('Supersaturation(nconc) (%)');
            ylabel('Height (km)');
            legend('C_ PI','C_ BG');
            title('SS from Fan calculation')
%             xlim([0,20]);
            ylim([0,15]);
            xlim([0,20]);
            hold off;
            
            subplot(2,2,4);
            hold on;
            line(ss_weighted.pi*100, ss_weighted.alt/1e3,'Linestyle','--', 'color','b','linewidth',2);
            line(ss_weighted.bg*100, ss_weighted.alt/1e3,'Linestyle','-', 'color','b','linewidth',2);
            xlabel('Weighted SS (%)');
            ylabel('Height (km)');
            legend('C_ PI','C_ BG');
%             xlim([0,15]);
            ylim([0,15]);
            xlim([0,20]);
            hold off;
            
        end
        
        function plot_slice_ss_latentheat()
            data = load('read_wrf_study_domain_fan.mat');
            height = data.height;
            lon = data.lon;
            lat = data.lat;
            altitude = zeros(size(height,3),1);
            PI = data.PI;
            BG = data.BG;
            %%%% find the index to filter the study domain
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            
          
            for i=1:size(height,3)
                height_cut = squeeze(height(:,:,i));
                altitude(i) = nanmean(height_cut(regionindx));
            end
            
            % slice at the 15th vertical layer, the height is 5km;
            layer = 43;
            time_indx = 30;
            ss_pi = nan(size(lon));
            ss_bg = nan(size(lat));
            
%             aa = nan(72,1);
%             for i=1:72
%                 aa(i) = max(PI.ss_fan(i,:,layer)*100);
%             end
            ss_pi(regionindx) = PI.ss(time_indx,:,layer)*100;
            ss_bg(regionindx) = BG.ss(time_indx,:,layer)*100;
            
            lc_pi = nan(size(lon));
            lc_bg = nan(size(lon));
            lc_pi(regionindx) = PI.lc_water(time_indx,:,layer);
            lc_bg(regionindx) = BG.lc_water(time_indx,:,layer);
            
            w_pi = nan(size(lon));
            w_bg = nan(size(lon));
            w_pi(regionindx) = PI.w(time_indx,:,layer);
            w_bg(regionindx) = BG.w(time_indx,:,layer);
            
            %%%filter out by w
%             ss_pi(ss_pi<0) = nan;
%             ss_bg(ss_bg<0) = nan;
%             
%             lc_pi(ss_pi<0) = nan;
%             lc_bg(ss_bg<0) = nan;
            
            figure;
            subplot(3,2,1);
            pcolor(lon,lat,ss_pi);
            xlim(lon_bdy);
            ylim(lat_bdy);
            colormap(blue_red_cmap);
            shading flat;
            colorbar;
            caxis manual;
            caxis([-10,10]);
            title('SS (%) in C PI');
            
            subplot(3,2,2);
            pcolor(lon,lat,ss_bg);
            xlim(lon_bdy);
            ylim(lat_bdy);
            colormap(blue_red_cmap);
            shading flat;
            caxis([-10,10]);
            colorbar;
            
            title('SS (%) in C BG');
            
            subplot(3,2,3);
            pcolor(lon,lat,lc_pi);
            xlim(lon_bdy);
            ylim(lat_bdy);
            colormap(blue_red_cmap);
            shading flat;
            colorbar;
            caxis([-2e-3,2e-3]);
            title('Latent heat rate (K/s) in C PI');
            
            subplot(3,2,4);
            pcolor(lon,lat,lc_bg);
            xlim(lon_bdy);
            ylim(lat_bdy);
            colormap(blue_red_cmap);
            shading flat;
            colorbar;
            caxis([-2e-3,2e-3]);
            title('Latent heat rate (K/s) in C BG');
            
            subplot(3,2,5);
            pcolor(lon,lat,w_pi);
            xlim(lon_bdy);
            ylim(lat_bdy);
            colormap(blue_red_cmap);
            shading flat;
            colorbar;
            caxis([-3,3]);
            title('Updraft velocity (m/s)');
            
            subplot(3,2,6);
            pcolor(lon,lat,w_bg);
            xlim(lon_bdy);
            ylim(lat_bdy);
            colormap(blue_red_cmap);
            shading flat;
            colorbar;
            caxis([-3,3]);
            title('Updraft velocity (m/s)');
            
        end
        
        function plot_wrf_pdf(var_name, var_bins)
            if strcmpi(var_name,'w')
                xlabelstr = 'Updraft velocity (m/s)';
            elseif strcmpi(var_name,'ss')
                xlabelstr = 'Supersaturation (%)';
            elseif strcmpi(var_name,'weighted_ss')
                xlabelstr = 'Supersaturation (%)';
            end
            
            
            
            data = load('read_wrf_study_domain.mat');
            height = data.height;
            lon = data.lon;
            lat = data.lat;
            altitude = zeros(size(height,3),1);
            PI = data.PI;
            BG = data.BG;
            %%%% find the index to filter the study domain
            lon_bdy = [-60.8, -60.5];
            lat_bdy = [-3.3, -3.1];
            regionindx = lon >= lon_bdy(1) & lon <= lon_bdy(2) & lat <= lat_bdy(2) & lat >= lat_bdy(1);
            
            for i=1:size(height,3)
                height_cut = squeeze(height(:,:,i));
                altitude(i) = nanmean(height_cut(regionindx));
            end
%             uv_bins = 2:0.1:24;
            al_bins = altitude/1e3;
            var_pdf.pi = nan(numel(var_bins), numel(altitude));
            var_pdf.bg = nan(numel(var_bins), numel(altitude));
            bin_width = (var_bins(2) - var_bins(1))/2;
            for i=1:numel(altitude)
                if strcmpi(var_name,'weighted_ss')
                    this_latent = squeeze(PI.latentrate(:,:,i));
                    this_var = squeeze(PI.ss(:,:,i)).*this_latent/sum(this_latent(:));
                else
                    this_var = squeeze(PI.(var_name)(:,:,i));
                end
                
                this_uv = squeeze(PI.w(:,:,i));
                this_var(this_uv<2) = nan;
                count = numel(this_var(:));
                if count > 0
                    for j=1:numel(var_bins)
                        indx = this_var <= var_bins(j)+bin_width & this_var > var_bins(j)-bin_width;
                        var_pdf.pi(j,i) = sum(indx(:))/count;
                    end
                end
                
                if strcmpi(var_name,'weighted_ss')
                    this_latent = squeeze(BG.latentrate(:,:,i));
                    this_var = squeeze(BG.ss(:,:,i)).*this_latent/sum(this_latent(:));
                else
                    this_var = squeeze(BG.(var_name)(:,:,i));
                end
                
                this_uv = squeeze(BG.w(:,:,i));
                this_var(this_uv<2) = nan;
                count = numel(this_var(:));
                if count > 0
                    for j=1:numel(var_bins)
                        indx = this_var <= var_bins(j)+bin_width & this_var > var_bins(j)-bin_width;
                        var_pdf.bg(j,i) = sum(indx(:))/count;
                    end
                end
                
            end
            
            var_pdf.pi(var_pdf.pi==0)=nan;
            var_pdf.bg(var_pdf.bg==0)=nan;
            
            [var_map,al_map] = meshgrid(var_bins*100, al_bins);
            figure;
            subplot(1,2,1);
            pcolor(var_map', al_map', var_pdf.pi*1e5);
            h = colorbar;
            caxis([0,65]);
            shading flat;
            ylim([0,15]);
            ylabel(h, 'Frequency of occurence (10^{-3} %)');
            ylabel('Height (km)');
            xlabel(xlabelstr);
            title('C_ PI');
            set(gca, 'fontsize', 14);
            
            subplot(1,2,2);
            pcolor(var_map', al_map', var_pdf.bg*1e5);
            h = colorbar;
            caxis([0,65]);
            shading flat;
            ylim([0,15]);
            ylabel('Height (km)');
            ylabel(h, 'Frequency of occurence (10^{-3} %)');
            xlabel(xlabelstr);
            title('C_ BG');
            set(gca, 'fontsize', 14);
            
        end
        
        function result = filter_fan(PI, BG, variable, altitude, opt)
            
            if strcmpi(variable,'weighted_ss')
                pi_var = PI.ss_fan;
                bg_var = BG.ss_fan;
            else
                pi_var = PI.(variable);
                bg_var = BG.(variable);
            end
            
            pi_w = PI.w;
            bg_w = BG.w;
            time = 14:1/12:20-1/12;
            
%             result.alt = [0:0.4:16]*1000;
            result.alt = altitude;
            switch opt
                case 'nofilter'
                    for i=1:numel(altitude)
                        this_pi_w = squeeze(pi_w(:,:,i));
                        indx = this_pi_w >2;
                        
                        this_pi_var = squeeze(pi_var(:,:,i));
                        if strcmpi(variable,'weighted_ss')
%                             this_pi_latent = PI.lc_water(:,:,i);
                            
                            this_pi_latent = PI.lc_water(:,:,i)+PI.lc_ice(:,:,i)+...
                                PI.lc_freeze(:,:,i)+PI.lc_rim(:,:,i);
%                             indx = indx&this_pi_latent>0;
                            result.pi(i) = nansum(this_pi_var(indx).*this_pi_latent(indx)...
                                /nansum(this_pi_latent(indx)));
                        else
                            result.pi(i) = nanmean(this_pi_var(indx));
                        end

                        this_bg_w = squeeze(bg_w(:,:,i));
                        this_bg_var = squeeze(bg_var(:,:,i));
                        indx = this_bg_w >2;
                        
                        if strcmpi(variable,'weighted_ss')
                            this_bg_latent = BG.lc_water(:,:,i)+BG.lc_ice(:,:,i)+...
                                BG.lc_freeze(:,:,i)+BG.lc_rim(:,:,i);
%                             this_bg_latent = BG.lc_water(:,:,i);
%                             indx = indx&this_bg_latent>0;
                            result.bg(i) = nansum(this_bg_var(indx).*this_bg_latent(indx)...
                                /nansum(this_bg_latent(indx)));
                        else
                            result.bg(i) = nanmean(this_bg_var(indx));
                        end

                    end
                    
                case 'alltime'
                    for i=1:numel(altitude)
                        this_pi_w = squeeze(pi_w(:,:,i));
                        this_pi_w(this_pi_w<2) = nan;
                        indx = this_pi_w > prctile(this_pi_w(:),90);
                        this_pi_var = squeeze(pi_var(:,:,i));
                        if strcmpi(variable,'weighted_ss')
%                             this_pi_latent = PI.lc_water(:,:,i)+PI.lc_ice(:,:,i)+...
%                                 PI.lc_freeze(:,:,i)+PI.lc_rim(:,:,i);
                            this_pi_latent = PI.lc_water(:,:,i);
                            result.pi(i) = nanmean(this_pi_var(indx).*this_pi_latent(indx)...
                                /nansum(this_pi_latent(indx)));
                        else
                            result.pi(i) = nanmean(this_pi_var(indx));
                        end

                        this_bg_w = squeeze(bg_w(:,:,i));
                        this_bg_w(this_bg_w<2) = nan;
                        this_bg_var = squeeze(bg_var(:,:,i));
                        indx = this_bg_w > prctile(this_bg_w(:), 90);
                        if strcmpi(variable,'weighted_ss')
%                             this_bg_latent = BG.lc_water(:,:,i)+BG.lc_ice(:,:,i)+...
%                                 BG.lc_freeze(:,:,i)+BG.lc_rim(:,:,i);
                            this_bg_latent = BG.lc_water(:,:,i);
                            result.bg(i) = nanmean(this_bg_var(indx).*this_bg_latent(indx)...
                                /nansum(this_bg_latent(indx)));
                        else
                            result.bg(i) = nanmean(this_bg_var(indx));
                        end

                    end
                    
                case 'avertime'
                    for i=1:numel(altitude)
                        this_pi_w = squeeze(pi_w(:,:,i));
                        this_pi_w(this_pi_w<2) = nan;
                        this_pi_var = squeeze(pi_var(:,:,i));
%                         if strcmpi(variable,'weighted_ss')
%                             this_pi_latent = nanmean(PI.lc_water(:,:,i),1);
%                             result.pi(i) = nansum(this_pi_var(this_pi_w > prctile(this_pi_w(:), 90 )).*this_pi_latent(this_pi_w > prctile(this_pi_w(:), 90 ))...
%                                 /nansum(this_pi_latent(this_pi_w > prctile(this_pi_w(:), 90 ))));
%                         else
%                             result.pi(i) = nanmean(this_pi_var(this_pi_w > prctile(this_pi_w(:), 90 )));
%                         end
                        result_time = nan(72,1);
                        for j=1:72
                            this_pi_var_time = this_pi_var(j,:);
                            this_pi_w_time = this_pi_w(j,:);
                            result_time(j) = nanmean(this_pi_var_time(this_pi_w_time>prctile(this_pi_w_time(:), 90 )));
                        end
                        result.pi(i) = nanmean(result_time);
                        
                        this_bg_w = squeeze(bg_w(:,:,i));
                        this_bg_w(this_bg_w<2) = nan;
                        this_bg_var = squeeze(bg_var(:,:,i));
%                         if strcmpi(variable,'weighted_ss')
%                             this_bg_latent = BG.lc_water(:,:,i);
%                             result.bg(i) = nansum(this_bg_var(this_bg_w > prctile(this_bg_w(:), 90 )).*this_bg_latent(this_bg_w > prctile(this_bg_w(:), 90 ))...
%                                 /nansum(this_bg_latent(this_bg_w > prctile(this_bg_w(:), 90 ))));
%                         else
%                             result.bg(i) = nanmean(this_bg_var(this_bg_w > prctile(this_bg_w(:), 90 )));
%                         end
                        result_time = nan(72,1);
                        for j=1:72
                            this_bg_var_time = this_bg_var(j,:);
                            this_bg_w_time = this_bg_w(j,:);
                            result_time(j) = nanmean(this_bg_var_time(this_bg_w_time>prctile(this_bg_w_time(:), 90 )));
                        end
                        result.bg(i) = nanmean(result_time);
                    end
                case 'averupdraft'
                    for i=1:numel(altitude)
                        temp_pi = [];
                        temp_bg = [];
                        for j=1:size(pi_w,1)
                            this_pi_w = squeeze(pi_w(j,:,i));
                                this_pi_w(this_pi_w<2) = nan;
                                this_pi_var = squeeze(pi_var(j,:,i));
                                if strcmpi(variable,'weighted_ss')
                                    this_pi_latent = PI.latentrate(j,:,i);
                                    temp_pi = [temp_pi; nanmean(this_pi_var(this_pi_w > prctile(this_pi_w(:), 90 )).*this_pi_latent(this_pi_w > prctile(this_pi_w(:), 90 ))...
                                        /nansum(this_pi_latent(this_pi_w > prctile(this_pi_w(:), 90 ))))];
                                else
                                    temp_pi = [temp_pi;nanmean(this_pi_var(this_pi_w > prctile(this_pi_w(:), 90 )))];
                                end

                                this_bg_w = squeeze(bg_w(j,:,i));
                                this_bg_w(this_bg_w<2) = nan;
                                this_bg_var = squeeze(bg_var(j,:,i));
                                if strcmpi(variable,'weighted_ss')
                                    this_bg_latent = BG.latentrate(j,:,i);
                                    temp_bg = [temp_bg; nanmean(this_bg_var(this_bg_w > prctile(this_bg_w(:), 90 )).*this_bg_latent(this_bg_w > prctile(this_bg_w(:), 90 ))...
                                        /nansum(this_bg_latent(this_bg_w > prctile(this_bg_w(:), 90 ))))];
                                else
                                    temp_bg = [temp_bg; nanmean(this_bg_var(this_bg_w > prctile(this_bg_w(:), 90 )))];
                                end
                        end
                        result.pi(i) = nanmean(temp_pi);
                        result.bg(i) = nanmean(temp_bg);
                    end
            end
            
            
%             result.pi = interp1(altitude, result.pi, result.alt);
%             result.bg = interp1(altitude, result.bg, result.alt);
%             
        end
        function bin_radius()
            a(1) = 2;
            for i=2:17
                a(i) = 2^(1/3)*a(i-1);
            end
        end
        
        function plot_v_ss_model()
            data = load('wrf_pi_14_bg');
            lon = data.lon;
            lat = data.lat;
            ss  = data.ss;
            vel = data.vel;
            tke = data.tke;
%             z = data.z;
            nc = data.nc;
            indx = lon(:,:,1) >-60.8 & lon(:,:,1) <-60.4 & lat(:,:,1) > -3.4 ...
                & lat(:,:,1) < -3.2;
            ss_l = [];
            vel_l = [];
            nc_l = [];
            tke_l = [];
%             z_l = [];
            for i=1:size(vel,3)-1
                for j=1:size(vel,4)
                    this_ss = squeeze(ss(:,:,i,j));
                    this_ss = this_ss(indx);
                    ss_l = [ss_l,this_ss(:)];
                    
                    this_vel = squeeze((vel(:,:,i,j)+vel(:,:,i+1,j))/2);
                    this_vel = this_vel(indx);
                    vel_l = [vel_l,this_vel(:)];
                    
                    this_nc = squeeze(nc(:,:,i,j));
                    this_nc = this_nc(indx);
                    nc_l = [nc_l,this_nc(:)];
                    
                    this_tke = squeeze(tke(:,:,i,j));
                    this_tke = this_tke(indx);
                    tke_l = [tke_l,this_tke(:)];
                end
            end
            
            figure;
%             indx = nc_l*1.225*1e-6>5;
%             ss_l = ss_l(indx);
%             vel_l = vel_l(indx);
%             nc_l = nc_l(indx);
            scatter(ss_l(:)*100,vel_l(:),[],nc_l(:)*1.225*1e-6,'filled');
            colormap('jet');
            h=colorbar;
            ylabel(h,'Nc [cm^{-3}]');
            xlabel('SS,%');
            ylabel('W, m/s');
            set(gca,'FontSize', 14);
%             xlim([-10,10]);
            line([-10,10],[0,0],'Color','black','LineStyle','--','Linewidth',1);
            line([0,0],[-10,10],'Color','black','LineStyle','--','Linewidth',1);
            
            figure;
            scatter(nc_l(:)*1.225*1e-6,tke_l(:));
        end
        % utility function
        function value = cut_by_region_indx(wrf_value,regionindx)
            wrf_value = permute(wrf_value, [4 3 1 2]);
            value = nan(size(wrf_value,1),sum(regionindx(:)), size(wrf_value,2));
            
            for i=1:size(wrf_value ,1)
                for j=1:size(wrf_value, 2)
                    this_layer = squeeze(wrf_value(i,j,:,:));
                    value(i,:,j) = this_layer(regionindx);
                end
                
            end
            
        end
    end
end