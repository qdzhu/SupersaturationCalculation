classdef misc_halo_analysis
    methods(Static = true)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property like methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = box_dir()
            value = '/Users/monicazhu/Box';
        end
        
        function value = halo_csv_dir()
            value = fullfile(misc_halo_analysis.box_dir, 'HALO_cleanup_csv');
        end
        
        function Dp = cdp_diameter()
            Dp = struct('low',[2.5, 2.9, 5, 7.5, 10.2, 11.8, 15.6, 18.7, 20.7, 24.6, 27.4, 29.2, 34.4, 39, 42.5],...
                        'up',[2.9, 5, 7.5, 10.2, 11.8, 15.6, 18.7, 20.7, 24.6, 27.4, 29.2, 34.4, 39, 42.5, 46],...
                        'mean',[2.69, 3.81, 6.12, 8.75, 10.97, 13.57, 17.08, 19.67, 22.57, 25.96, 28.29, 31.69, 36.63, 40.71, 44.22]);
        end
        
        function Dp = cas_diameter()
            % remaining issue, the first two bins may be aerosols instead
            % of cloud droplet
            Dp = struct('low',[0.89,0.96,3,5,7.2,15,20,25,30,35,40,45],...
                        'up',[0.96,3,5,7.2,15,20,25,30,35,40,45,50]);
                    
        end
        %%%%%%%%%%%%%%%%%%%
        % Utility methods %
        %%%%%%%%%%%%%%%%%%%
        
        % return CDPFILE structure including two fields, filename
        % and the date extracted from the filename
        function cdpfile = read_cdp_file()
            filepattern = '*CCP_CDP*';
            cdpdir = dir(fullfile(misc_halo_analysis.halo_csv_dir, filepattern));
            cdpfile = make_empty_struct_from_cell({'name','date'});
            for i=1:numel(cdpdir)
                cdpfile(i).('name') = fullfile(misc_halo_analysis.halo_csv_dir, cdpdir(i).name);
                cdpfile(i).('date') = datenum(misc_halo_analysis.date_from_halo_filename(cdpdir(i).name));
            end
        end
        
        function casfile = read_cas_file()
            filepattern = '*CAS_DPOL*';
            casdir = dir(fullfile(misc_halo_analysis.halo_csv_dir, filepattern));
            casfile = make_empty_struct_from_cell({'name','date'});
            for i=1:numel(casdir)
                casfile(i).('name') = fullfile(misc_halo_analysis.halo_csv_dir, casdir(i).name);
                % time formate dd_mm_yyyy
                thisdate = regexp(casdir(i).name, '\d\d_\d\d_\d\d\d\d','match','once');
                casfile(i).('date') = datenum(strcat(thisdate(7:10),'-',thisdate(4:5),'-',thisdate(1:2)));
            end
        end
        
        function adlrfile = read_adlr_file()
            filepattern = '*adlr*';
            adlrdir = dir(fullfile(misc_halo_analysis.halo_csv_dir, filepattern));
            adlrfile = make_empty_struct_from_cell({'name','date'});
            for i=1:numel(adlrdir)
                adlrfile(i).('name') = fullfile(misc_halo_analysis.halo_csv_dir, adlrdir(i).name);
                adlrfile(i).('date') = datenum(misc_halo_analysis.date_from_halo_filename(adlrdir(i).name));
            end
        end
        
        function date = date_from_halo_filename(filename)
            % the date formate is yyyymmdd
            dstr = regexp(filename, '\d\d\d\d\d\d\d\d','match','once');
            date = strcat(dstr(1:4),'-',dstr(5:6),'-',dstr(7:8));
        end
        
        function matchfile = read_file_match()
            cdpfile = misc_halo_analysis.read_cdp_file;
            casfile = misc_halo_analysis.read_cas_file;
            adlrfile = misc_halo_analysis.read_adlr_file;
            
            matchfile = make_empty_struct_from_cell({'cdpname','casname','adlrname','date'});
            n = numel(cdpfile); %% cdp has the least observation days
            
            for i=1:n
                matchfile(i).('cdpname') = cdpfile(i).('name');
                matchfile(i).('date') = datestr(cdpfile(i).('date'));
                matchfile(i).('casname') = casfile(extractfield(casfile,'date') == cdpfile(i).('date')).('name');
                matchfile(i).('adlrname') = adlrfile(extractfield(adlrfile,'date') == cdpfile(i).('date')).('name');
            end
            
            
        end
        
        
        
        function time = common_utcsec(t1,t2,t3)
            time = unique([t1,t2,t3]);
            indx = false(size(time));
            for i=1:numel(time)
                if ismember(time(i), t1) && ismember(time(i), t2) && ismember(time(i), t3)
                    indx(i) = true;
                end
            end
            time = time(indx);
        end
        
        function cdpdata = read_cdp_data(cpdname)
            file = csvread(cpdname,3,0);
            cdpdata = make_empty_struct_from_cell({'utcsec','meandp','nconc'});
            cdpdata.utcsec = fix(file(:,1));
            cdpdata.meandp = file(:,2)/2;% convert from diameter to radius
            
            lowdp = misc_halo_analysis.cdp_diameter.low;
            updp = misc_halo_analysis.cdp_diameter.up;
            dlogdp = log(updp)-log(lowdp);
            nconc = zeros(size(cdpdata.meandp));
            for i_bin = 1:numel(lowdp)
                bin_indx = i_bin+3;
                nconc = nconc + file(:,bin_indx)./file(:,3); %*dlogdp(i_bin)./file(:,3);
            end
            cdpdata.nconc = nconc;
            % remove the filling value
            cdpdata.meandp(cdpdata.meandp >10000 | cdpdata.meandp <0 )= nan;
            cdpdata.nconc(cdpdata.nconc > 10000 | cdpdata.nconc < 0) = nan;
        end
        
        function casdata = read_cas_data(casname)
            file = csvread(casname, 3, 0);
            casdata = make_empty_struct_from_cell({'utcsec','meandp','nconc','lwc'});
            casdata.utcsec = fix(file(:,1));
            casdata.meandp = zeros(size(casdata.utcsec));
            casdata.nconc = zeros(size(casdata.utcsec));
            casdata.lwc = file(:,end-3);
            casdata.lwc(casdata.lwc<0) = nan;
            lowdp = misc_halo_analysis.cas_diameter.low;
            updp = misc_halo_analysis.cas_diameter.up;
            meandp = (lowdp+updp)/2/2;%convert from diameter to radius
            nconc = zeros(size(casdata.utcsec));
            meandp_time_nconc = zeros(size(casdata.utcsec));
            for i_bin = 3:numel(lowdp)
                bin_indx = i_bin+1;
                nconc = nconc+file(:,bin_indx);
                meandp_time_nconc = meandp_time_nconc + file(:,bin_indx)*meandp(i_bin);
                % ignore the missing value
                indx = file(:,bin_indx) < 0;
                nconc(indx) = nan;
                meandp_time_nconc(indx) = nan;
            end
            casdata.meandp = meandp_time_nconc./nconc;
            casdata.nconc = nconc;
        end
        
        function adlrdata = read_adlr_data(adlrname)
            file = csvread(adlrname, 3, 0);
            adlrdata = make_empty_struct_from_cell({'utcsec','temp','w','alt'});
            adlrdata.utcsec = fix(file(:,1));
            adlrdata.alt = file(:,5);
            adlrdata.temp = file(:,21);
            adlrdata.w = file(:,18);
            adlrdata.w(adlrdata.w<=-100) = nan;
            adlrdata.temp(adlrdata.temp<=-100) = nan;
        end
        
        function struc_field = struct_filter(istruc, target_sec, fieldname)
            istruc_field = extractfield(istruc, fieldname);
            istruc_utcsec = extractfield(istruc, 'utcsec');
            struc_field = nan(size(target_sec));
            for i=1:numel(target_sec)
                indx = istruc_utcsec == target_sec(i);
                struc_field(i) = nanmean(istruc_field(indx));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Make methods    %
        %%%%%%%%%%%%%%%%%%%
        
        % return a data struct including all necessary fields for future
        % calculation/analysis
        function make_match_data()
            matchfile = misc_halo_analysis.read_file_match;
            match = make_empty_struct_from_cell({'utcsec','alt','w','temp','lwc','cas_meandp','cas_nconc','cdp_meandp','cdp_nconc','date'});
            for i_date = 1:numel(matchfile)
                cdpdata = misc_halo_analysis.read_cdp_data(matchfile(i_date).cdpname);
                casdata = misc_halo_analysis.read_cas_data(matchfile(i_date).casname);
                adlrdata = misc_halo_analysis.read_adlr_data(matchfile(i_date).adlrname);
                % do data filter based on the utcsec
                match(i_date).date = matchfile(i_date).date;
                match(i_date).utcsec = misc_halo_analysis.common_utcsec(extractfield(cdpdata,'utcsec'), extractfield(casdata,'utcsec'), extractfield(adlrdata,'utcsec'));
                
                %cdp_indx = ismember(cdpdata.utcsec, match(i_date).utcsec);
                match(i_date).cdp_meandp = misc_halo_analysis.struct_filter(cdpdata, match(i_date).utcsec, 'meandp');
                match(i_date).cdp_nconc = misc_halo_analysis.struct_filter(cdpdata, match(i_date).utcsec, 'nconc');
                %cas_indx = ismember(casdata.utcsec, match(i_date).utcsec);
                match(i_date).cas_meandp = misc_halo_analysis.struct_filter(casdata, match(i_date).utcsec, 'meandp');
                match(i_date).cas_nconc = misc_halo_analysis.struct_filter(casdata, match(i_date).utcsec,'nconc');
                match(i_date).lwc = misc_halo_analysis.struct_filter(casdata, match(i_date).utcsec, 'lwc');
                %adlr_indx = ismember(adlrdata.utcsec, match(i_date).utcsec);
                match(i_date).temp = misc_halo_analysis.struct_filter(adlrdata, match(i_date).utcsec,'temp');
                match(i_date).w = misc_halo_analysis.struct_filter(adlrdata, match(i_date).utcsec, 'w');
                match(i_date).alt = misc_halo_analysis.struct_filter(adlrdata, match(i_date).utcsec,'alt');
            end
            save('halo_match.mat','match');
        end
        
        function make_supersaturation()
            data = load('halo_match.mat');
            match = data.match;
            [a0,a1,a2,a3,a4,a5,a6,Rg,Ra,Cpa,Mma,Rv,Cpv,Mmv,pl,ps,Mms,alpha,w,Po,To,g,k_mu,k_ml]=recalculate_CAIPEEX_result.Constant;
            for i_date=1:numel(match)
                 T = match(i_date).temp;
                 vel = match(i_date).w;
                 Tc = T-To;
                 Ho = match(i_date).alt;
                 P = Po*exp(-g*Ho./(Ra*T));
                 L=2.495e6-2.3e3*Tc;          % latent heat of evaporation
                 D=(2.26e-5+1.5e-7*Tc)*Po./P;  % diffusion coeff.
                 A = (g*L./(Cpa*Rv*T.^2)-g./(Ra*T));
                 cdp_meandp = match(i_date).cdp_meandp;
                 cdp_nconc = match(i_date).cdp_nconc;
                 cdp_meandp(cdp_meandp==0) = nan; % avoid inf SS
                 cdp_nconc(cdp_nconc==0) = nan;
                 match(i_date).cdp_SS = A.*vel./(4*pi*D.*cdp_meandp.*cdp_nconc)*100;
                 cas_meandp = match(i_date).cas_meandp;
                 cas_nconc = match(i_date).cas_nconc;
                 cas_meandp(cas_meandp == 0) =nan;
                 cas_nconc(cas_nconc ==0) = nan;
                 match(i_date).cas_SS = A.*vel./(4*pi*D.*cas_meandp.*cas_nconc)*100;
            end
            save('halo_match_ss.mat','match');
        end
        
        %%%%%%%%%%%%%%%%%%%
        % Plot methods    %
        %%%%%%%%%%%%%%%%%%%
        function plot_updraft_wind_supersaturation_oneday(n)
            data = load('halo_match_ss.mat');
            match = data.match(n);
            
            figure;
            indx = match.cdp_nconc>1;
            scatter(match.cdp_SS(indx), match.w(indx),[],match.lwc(indx), 'filled');
            h = colorbar;
            caxis([0,0.1])
            xlabel('SS,%');
            ylabel('W, ms^{-1}');
            ylabel(h,'LWC, g cm^{-3}');
            
            figure;
            hold on;
            subplot(2,2,1);
            line(match.cdp_meandp*2, match.cas_meandp*2,'linestyle','none','marker','o','markersize',2,'color','blue');
            line([0,50],[0,50],'linestyle','--','color','red','linewidth',2);
            xlabel('CDP');
            ylabel('CAS');
            title('Mean Diameter (microns)');
            
            subplot(2,2,2);
            line(match.cdp_nconc, match.cas_nconc,'linestyle','none','marker','o','markersize',2,'color','blue');
            line([0,2500],[0,2500],'linestyle','--','color','red','linewidth',2);
            xlabel('CDP');
            ylabel('CAS');
            title('Number Concentration (#cm^{-3})');
            
            subplot(2,2,3);
            line(match.cdp_SS, match.cas_SS,'linestyle','none','marker','o','markersize',2,'color','blue');
            line([0,0],[-100,100],'linestyle','--','color','red','linewidth',2);
            line([-100,100],[0,0],'linestyle','--','color','red','linewidth',2);
            xlabel('CDP');
            ylabel('CAS');
            xlim([-100,100]);
            ylim([-100,100]);
            title('Supersaturation (%)');
            
            
            
            figure;
            %match.cdp_nconc(match.cdp_nconc>10) = nan;
            [~,edges] = histcounts(log10(match.cdp_nconc(indx)));
             histogram(match.cdp_nconc(indx),10.^edges)
            set(gca, 'xscale','log')
            xlabel('CDP Nconc');
            ylabel('frequencey');
            
            figure;
            histogram(match.cdp_meandp(indx)*2);
            xlabel('CDP Mean diameter');
            ylabel('frequencey');
            
            figure;
            [~,edges] = histcounts(log10(match.lwc(indx)));
            histogram(match.lwc(indx),10.^edges)
            set(gca, 'xscale','log')
            xlabel('CDP LWC (g cm^{-3})');
            ylabel('frequencey');
            
            figure;
            cloud_edge_indx = match.lwc(indx) < prctile(match.lwc(indx),5);
            cloud_non_edge_indx = match.lwc(indx)>= prctile(match.lwc(indx),5);
            cdp_SS = match.cdp_SS(indx);
            ss_edge = cdp_SS(cloud_edge_indx);
            ss_nonedge = cdp_SS(cloud_non_edge_indx);
            subplot(1,2,1);
            hold on;
            histogram(ss_edge);
            hold off;
            
            subplot(1,2,2);
            hold on;
            histogram(ss_nonedge);
            hold off;
        end
        
        function plot_hist_side_plot(x,y)
            trace1 = struct(...
                'x', x, ...
                'y', y, ...
                'mode', 'markers', ...
                'name', 'points', ...
                'marker', struct(...
                'color', 'rgb(102,0,0)', ...
                'size', 2, ...
                'opacity', 0.4), ...
                'type', 'scatter');
            trace2 = struct(...
                  'x', x, ...
                  'name', 'x density', ...
                  'marker', struct('color', 'rgb(102,0,0)'), ...
                  'yaxis', 'y2', ...
                  'type', 'histogram');
            trace3 = struct(...
                  'y', y, ...
                  'name', 'y density', ...
                  'marker', struct('color', 'rgb(102,0,0)'), ...
                  'xaxis', 'x2', ...
                  'type', 'histogram');  
            data = {trace1, trace2, trace3};
            layout = struct(...
                'showlegend', false, ...
                'autosize', false, ...
                'width', 600, ...
                'height', 550, ...
                'xaxis', struct(...
                  'domain', [0, 0.85], ...
                  'showgrid', false, ...
                  'zeroline', false), ...
                'yaxis', struct(...
                  'domain', [0, 0.85], ...
                  'showgrid', false, ...
                  'zeroline', false), ...
                'margin', struct('t', 50), ...
                'hovermode', 'closest', ...
                'bargap', 0, ...
                'xaxis2', struct(...
                  'domain', [0.85, 1], ...
                  'showgrid', false, ...
                  'zeroline', false), ...
                'yaxis2', struct(...
                  'domain', [0.85, 1], ...
                  'showgrid', false, ...
                  'zeroline', false));
             response = plotly(data, struct('layout', layout, 'fileopt', 'overwrite'));
        end
    end
end