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
        
        % return a data struct including all necessary fields for future
        % calculation/analysis
        function data = read_match_data()
            matchfile = misc_halo_analysis.read_file_match;
            for i_date = 1:numel(matchfile)
                cdpdata = misc_halo_analysis.read_cdp_data(matchfile(i_date).cdpname);
                casdata = misc_halo_analysis.read_cas_data(matchfile(i_date).casname);
                adlrdata = misc_halo_analysis.read_adlr_data(matchfile(i_date).adlrname);
            end
            
            % do data filter based on the utcsec
        end
        
        function cdpdata = read_cdp_data(cpdname)
            file = csvread(cpdname,3,0);
            cdpdata = make_empty_struct_from_cell({'utcsec','meandp','nconc'});
            cdpdata.utcsec = file(:,1);
            cdpdata.meandp = file(:,2);
            
            lowdp = misc_halo_analysis.cdp_diameter.low;
            updp = misc_halo_analysis.cdp_diameter.up;
            dlogdp = log(updp)-log(lowdp);
            nconc = zeros(size(cdpdata.meandp));
            for i_bin = 1:numel(lowdp)
                bin_indx = i_bin+3;
                nconc = nconc + file(:,bin_indx)*dlogdp(i_bin);
            end
            cdpdata.nconc = nconc;
        end
        
        function casdata = read_cas_data(casname)
            file = csvread(casname, 3, 0);
            casdata = make_empty_struct_from_cell({'utcsec','meandp','nconc','lwc'});
            casdata.utcsec = file(:,1);
            casdata.meandp = zeros(size(casdata.utcsec));
            casdata.nconc = zeros(size(casdata.utcsec));
            casdata.lwc = file(:,end-3);
            casdata.lwc(casdata.lwc<0) = nan;
            lowdp = misc_halo_analysis.cas_diameter.low;
            updp = misc_halo_analysis.cas_diameter.up;
            meandp = (lowdp+updp)/2;
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
            adlrdata.utcsec = file(:,1);
            adlrdata.alt = file(:,5);
            adlrdata.temp = file(:,12);
            adlrdata.w = file(:,18);
            adlrdata.w(adlrdata.w<=-100) = nan;
            adlrdata.temp(adlrdata.temp<=-100) = nan;
        end
        %%%%%%%%%%%%%%%%%%%
    end
end