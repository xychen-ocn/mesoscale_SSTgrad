% This script will look for "feature based" SST gradient that is associated
% with mesoscale varibility in the skin temperature derived from the
% RHB seasnake measurement. 
% mesoscale: 50km - 100km scale
% this scale will be revealed via the scalogram produced from wavelet
% analysis. 


%% input data: rhb_group . (data grouped by day);

%% data processing steps:
%  - 1. 10min SST samples --> interpolate to be a 2km uniformally spaced
%  dataset. 
%  - 2. detrend and low pass filter --> filter out varaibility smaller than 20km;
%  - 3. find peaks and the length scale of the variation associated with
%  at each sample location. 
%  - 4. manually segmentize the SST records to pick out realively strong
%  "frontal" events (dT>0.25°C)
%

%% saving data:
%  - 1. save data by day and segment number. 
clear all; clc; close all;

%% ------------------------- Code Part -------------------------------- %%
loc_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
addpath([loc_datadir filesep '../bin']);

% --- a. load in data ---- %
RHB_dataFN = 'rhb_daily_grouped_10min_data_0909latest.mat';
load([loc_datadir filesep RHB_dataFN]);

% --- b. choose to work with the 8 Strong SST variance days --- %
added_days = [datenum('20200114','yyyymmdd'), datenum('20200130','yyyymmdd')];
DOIs = [rhbdates.strong_SSTvar added_days];                                             % Date of Interests.


% --- c. for each day, execute the data processing steps and segmentize the
% ---    SST gradient:
spatial_res = 2;                                                           % units: km/sample  
wvlen_cutoff =20;                                                          % units: km, for filtering
VOIs = {'local_time','time','lon','lat','sst_seasnake', 'latent_heat_flux', ...
       'sensible_heat_flux','tau','hnet','wspd_10N','wdir','hl_bulk','hs_bulk', ...
       'MO_length','wave_flag','ship_contamination_flag','QC_flag'};

%DOIs = rhbdates.moderate_SSTvar;

for id =1:length(DOIs)
    DN = datestr(DOIs(id), 'mmmdd');
    if ~exist([loc_datadir filesep DN] ,'dir')
        mkdir(loc_datadir, DN);
    end
    %    seg = struct([]);                                                          % initialize the segment structure.
    
    
    %%%%%%%%--> step 1:
    moving_flag = rhb_group.(DN).in_transit_flag;
    
    if id ==8
        % make the end of the record moving because we know it is
        % transiting a front.
        ids_tmp =find(moving_flag);
        
         moving_flag(max(ids_tmp):end)=true;
    end
    
    sst = rhb_group.(DN).sst_seasnake(moving_flag);
    wspd_10N =rhb_group.(DN).wspd_10N(moving_flag);
    traj = rhb_group.(DN).traj(moving_flag);
    dist_equal = [0:spatial_res:round(max(traj))];                         % units: km
    SST_interp = interp1(traj, sst ,dist_equal,'linear', 'extrap');
    U10N_interp = interp1(traj, wspd_10N, dist_equal, 'linear','extrap');
    
    %%%%%%%%--> step 2:
    %%% -- detrend:
    SST_interp_dtr = detrend(SST_interp,1);                                % detrend (I don't need this.)
    trend_linear = SST_interp - SST_interp_dtr;
    
    U10N_interp_dtr = detrend(U10N_interp,1);
    trend_linear_U10N = U10N_interp - U10N_interp_dtr;
    
    %%% -- low-pass filter:
    % wpass: normalized passband frequency: pi rad/sample.
    % filter out wavelength smaller than 20km;
    % wavenumber: 2*pi rad/20km;  resolution: 2km per sample
    % rad per sample = wavenumber_cutoff * resolution
    wn_cutoff = 2/wvlen_cutoff ;  % pi rad/km
    wpass = wn_cutoff * spatial_res;
    SST_lowpassed_dtr = lowpass(SST_interp_dtr, wpass);   % retain frequency lower than the cutoff;
    SST_lowpassed = SST_lowpassed_dtr+trend_linear;
    
    U10N_lowpassed_dtr = lowpass(U10N_interp_dtr, wpass);
    U10N_lowpassed = U10N_lowpassed_dtr + trend_linear_U10N;
    
    %%%%%%%%--> step 3:
    if id==10
        minpeakdist = 18;
    else
        minpeakdist = 20;
    end
    [pks.val, pks.loc, pks.w, pks.p] = findpeaks(SST_lowpassed, dist_equal, 'MinPeakDistance',minpeakdist, 'MinPeakProminence',0.05);
    [trghs.val, trghs.loc, trghs.w, trghs.p] = findpeaks(-(SST_lowpassed), dist_equal, 'MinPeakDistance',20,'MinPeakProminence',0.01);
    [wt,fn]=cwt(SST_interp,'morse');
    %[wt,fn]=cwt(SST_interp,'bump');
    % [wt,fn]=cwt(SST_interp,'amor');
    
    % -- produce a sanity check figure:
    mag_wt = abs(wt);
    wvlen = 1./fn .* 2;    % 1cycle/ (cycle per sample) * spatial resolution (km/sample) = km
    
    [DD, WL] = meshgrid(dist_equal, wvlen);
    BW = imregionalmax(abs(wt));
    % take the maximum that has a wavelength > 50km
    C1 = WL>25;
    C2 = mag_wt>0.75*max(abs(wt(:)));     % magnitude criterion set to 75% of the maximum magnitude.
    cond = (C1 & C2 & BW);
    
    % lengthscale of the local maximum :
    wvlen_locmax = WL(cond);
    xloc_locmax = DD(cond);
    
    %% plot
    hfig = figure(10);clf;
    %%========== signal ========= %
    subplot(2,2,1)
    %yyaxis left
    plot(dist_equal,SST_interp,'linewidth',1.1)
    hold on
    plot(traj, sst,'.b');
    plot(dist_equal,trend_linear,'--k');
    xlabel('km');
    ylabel('Celcius');
    set(gca,'fontsize',14);
    
    %yyaxis right
    %plot(dist_equal,detrend(SST_interp,0),'o');
    hold on
    plot(pks.loc, pks.val,'*m','linewidth',1.2, 'markersize',12);
    plot(trghs.loc, -trghs.val,'*b','linewidth',1.2, 'markersize',12);
    plot(dist_equal, SST_lowpassed,'-r','linewidth',1.2);
    %     ylabel('Celcius')
    %     set(gca,'fontsize',14);
    xrange_p1 = get(gca,'xlim');
    yrange = get(gca,'ylim');
    ylim([yrange(1), yrange(1)+1.2]);
    grid on
    title(DN)
    %hold off
    
    %%====== Scalogram ======== %%
    subplot(2,2,3)
    levs = max(mag_wt(:)).*[0.5:0.1:1];
    pcolor(dist_equal, wvlen , mag_wt);
    shading flat
    hold on
    plot(DD(cond), WL(cond),'om','MarkerSize',10);
    [c,h] = contour(dist_equal, wvlen, mag_wt, levs,'k');
    axis tight
    xlabel('km');
    ylabel('wavelengths (km)');
    set(gca,'yscale','linear');
    hb = colorbar;
    set(hb,'location','southout','orientation','horizontal');
    set(get(hb,'xlabel'),'String','magnitude');
    set(gca,'ytick',[5, 10, 50,100],'TickDir','both')
    grid on
    xlim(xrange_p1 )
    set(gca,'fontsize',14);
    title('Magnitude Scalogram');
    hold off
    
    
    %%%%%%%%--> step 4:
    % compute gradient using the selected peaks and the consecutive
    % troughs:
    %
    % -- add: focus on pks that have a magnitude in the wavelet
    % transformation > 0.75 % of the maximum magnitude:
    
    
    seg = find_fronts_and_compute_gradient(rhb_group.(DN), moving_flag, dist_equal, SST_lowpassed, pks, trghs, mag_wt, WL,traj, VOIs, hfig);
    
    
    figure(10);
    subplot(2,2,[2,4]);
    hold on;
    plot(rhb_group.(DN).lon, rhb_group.(DN).lat,'--b');
    plot(rhb_group.(DN).lon(1), rhb_group.(DN).lat(1),'*k','markersize',13);
    plot(rhb_group.(DN).lon(end), rhb_group.(DN).lat(end),'dk','markersize',13);
    
    % overlay the surface wind direction:
    wndir = get_cartesian_direction(rhb_group.(DN).wdir, 'Meteo');
    u = rhb_group.(DN).wspd_10N(moving_flag).* cosd(wndir(moving_flag));
    v = rhb_group.(DN).wspd_10N(moving_flag).* sind(wndir(moving_flag));
    %scatter(rhb_group.(DN).lon, rhb_group.(DN).lat, 10, rhb_group.(DN).local_time);
    
    if median(wndir)>180
        scale =-1;
    else
        scale=1;
    end
    scatter(rhb_group.(DN).lon(moving_flag), rhb_group.(DN).lat(moving_flag)+ scale*0.025, 20, rhb_group.(DN).wspd_10N(moving_flag),'filled');
    hb2=colorbar;
    set(get(hb2,'xlabel'),'string','m/s');
    %datetick(hb2,'x','mmmdd-hhZ','keepticks');
    caxis([median(rhb_group.(DN).wspd_10N)-2*std(rhb_group.(DN).wspd_10N), ...
        median(rhb_group.(DN).wspd_10N)+2*std(rhb_group.(DN).wspd_10N)])
    %     quiver(rhb_group.(DN).lon(1:3:end), rhb_group.(DN).lat(1:3:end)-0.025, ...
    %         u(1:3:end), v(1:3:end), 0.5,'k');
    rhb_lon = rhb_group.(DN).lon(moving_flag);
    rhb_lat = rhb_group.(DN).lat(moving_flag);
    
    quiver(rhb_lon(1:3:end), rhb_lat(1:3:end)+ scale*0.025, ...
        u(1:3:end), v(1:3:end), 0.5,'k');
    axis('square')
    xlabel('longitude');
    ylabel('latitude');
    set(gca,'fontsize',14)
    
    if 1==0
        cnt = 0;
        num_pks = length(pks.val);
        for ip = 1:num_pks
            
            % get the mag_wt at the xlocs of the pks:
            xid = find(dist_equal==pks.loc(ip));
            TF = islocalmax(mag_wt(:,xid));
            [feature_wvlen, maxid] = max(WL(TF,xid));
            good_peaks = (max(mag_wt(TF,xid))>0.70*max(mag_wt(:))) & (round(feature_wvlen)>=30);
            
            if good_peaks || pks.p(ip)> min(pks.p)
                % enter manual selection:
                % show current peak from figure 10.
                figure(10)
                subplot(2,2,1)
                plot(pks.loc(ip), pks.val(ip),'^m','linewidth',1.5, 'markersize',13.5);
                
                
                % manually select.
                prompt = 'compute gradient (Y/N)?';
                Ystr = input(prompt,'s');
                
                
                if strcmpi(Ystr,'Y')
                    if cnt ==0
                        n = 1;
                    else
                        n = cnt+1;
                    end
                    
                    % compute gradient:
                    % ---  i. find the consecutive troughs:
                    dx = abs(pks.loc(ip) - trghs.loc);
                    
                    % take the first two.
                    [dx_ascending, sid]= sort(dx);
                    qttys = fieldnames(trghs);
                    for iv =1 :4
                        tmp = trghs.(qttys{iv})(sid);
                        trghs_sel.(qttys{iv}) = tmp(1:2);
                    end
                    
                    % remove troughs that does not satisfy the following condition:
                    % require dx to be larger than 12.5km.
                    dx_upperlimit = dx_ascending(1:2)<=feature_wvlen;
                    dx_lowerlimit = dx_ascending(1:2)>12.5;
                    cond_dist = dx_lowerlimit & dx_upperlimit;
                    for iv =1:4
                        trghs_sel.(qttys{iv}) = trghs_sel.(qttys{iv})(cond_dist);
                    end
                    
                    
                    % if this is the last peak & only 1 trough is selected
                    if (ip == num_pks) & (length(trghs_sel.loc)==1)
                        % add the end point in the SST record to compute
                        % gradient.
                        % if there is no trough to select:
                        if pks.loc(ip)>trghs.loc(end)
                            if ((pks.val(ip)-SST_lowpassed(end))>0.1)
                                if ((dist_equal(end) - pks.loc(ip))> 0.25*feature_wvlen) &  ...
                                        ((dist_equal(end) - pks.loc(ip))<=0.8*feature_wvlen)
                                    trghs_sel.loc(2) = dist_equal(end);
                                    trghs_sel.val(2) = -SST_lowpassed(end);
                                else
                                    trghs_sel.loc(2) = pks.loc(ip)+ 0.5*feature_wvlen;
                                    trghs_sel.val(2) = - interp1(dist_equal,SST_lowpassed, trghs_sel.loc(2),'linear');
                                end
                            end
                        else
                            % if there is a trough to select but exceed upperlimit:
                            if find(dx_upperlimit)==1
                                % find out on which side of the peak:
                                sideID = find(~dx_upperlimit);
                                trghs_sel.loc(2) = pks.loc(ip)+0.5*feature_wvlen;
                                trghs_sel.val(2) = -interp1(dist_equal, SST_lowpassed, trghs_sel.loc(2),'linear');
                                
                                % if this make rest of the section long enough
                                % to compute gradient between the added second
                                % point and the troughs: flag it and compute
                                % the gradient later on:
                                dist_remained = dx_ascending(sideID) - 0.5*feature_wvlen;
                                trwvlen = find_feature_wavelen(trghs.loc(end), dist_equal, mag_wt, WL);
                                
                                if dist_remained <trwvlen
                                    flag_revisit_last_trgh=true;
                                end
                            end
                        end
                    end
                    
                    % if this is the first peak  & only 1 trough is detected:
                    if (ip==1) & (length(trghs_sel.loc)==1)
                        % add the first point in the SST record to compute
                        % gradient:
                        if ((pks.val(ip)-SST_lowpassed(1))>0.1)
                            if ((pks.loc(ip)-dist_equal(1))> 0.25*feature_wvlen) &  ...
                                    ((pks.loc(ip)-dist_equal(1))<=0.8*feature_wvlen)
                                trghs_sel.loc(2) = dist_equal(1);
                                trghs_sel.val(2) = -SST_lowpassed(1);
                            elseif  ((pks.loc(ip)-dist_equal(1))>0.8*feature_wvlen)
                                trghs_sel.loc(2) = pks.loc(ip)- 0.5*feature_wvlen;
                                trghs_sel.val(2) = - interp1(dist_equal,SST_lowpassed, trghs_sel.loc(2),'linear');
                                
                            end
                        end
                    end
                    
                    % if ip is not the first or the last but still only 1
                    % trough is selected:
                    if (ip>1&ip<num_pks) & (length(trghs_sel.loc)<2)
                        % if this is due to an upper limit issue:
                        if find(dx_upperlimit)==1
                            % find out on which side of the peak:
                            sideID = find(~dx_upperlimit);
                            trghs_sel.loc(2) = pks.loc(ip)+0.5*feature_wvlen;
                            trghs_sel.val(2) = -interp1(dist_equal, SST_lowpassed, trghs_sel.loc(2),'linear');
                        end
                    end
                    
                    
                    % sanity check:
                    figure(10)
                    subplot(2,2,1)
                    hold on
                    plot(trghs_sel.loc, -trghs_sel.val,'ob','linewidth',1.5, 'markersize',13.5);
                    
                    sst_dif = pks.val(ip) - (-trghs_sel.val);
                    sst_grad_abs = sst_dif ./ (pks.loc(ip) - (trghs_sel.loc)) ;
                    
                    
                    xlocs = [trghs_sel.loc];
                    [xlocs_sorted, sxid] = sort(xlocs);
                    sst_grad_abs = sst_grad_abs(sxid);
                    sst_dif = sst_dif(sxid);
                    
                    % figure out how the ship is travel relative to the wind
                    % from trough to peak and from peak to troughs.
                    % I can find this out from a mask:
                    ship_wind_align = rhb_group.(DN).ship_wind_align_mask_interp;
                    
                    % find the wavelength revealed from the scalogram:
                    %                 xid = find(dist_equal ==pks.loc(ip));
                    %                 TF = islocalmax(mag_wt(:,xid));
                    %                 feature_wvlen = round(max(WL(TF,xid)));
                    
                    %% if gradient is taken from two sides of the peak
                    if length(xlocs_sorted)==2
                        % segment 01: trough 1 -> peak0
                        xmask = (traj>=xlocs_sorted(1)) & (traj<=pks.loc(ip));          % traj has the same length size as the time-based variables.
                        for iv = 1:length(VOIs)
                            varn = VOIs{iv};
                            seg(n).(varn) = rhb_group.(DN).(varn)(xmask);
                        end
                        xmask2 = (dist_equal>=xlocs_sorted(1)) & (dist_equal<=pks.loc(ip));
                        seg(n).SSTgrad = sst_grad_abs(1) .* mode(sign(ship_wind_align(xmask2)));
                        seg(n).front_loc_traj = 0.5*(xlocs_sorted(1)+ pks.loc(ip));        % km;
                        % translate back to the lat-lon coordinate.
                        [minval, minid] = min(abs(traj - 0.5*(xlocs_sorted(1)+ pks.loc(ip))));
                        front_lon = rhb_group.(DN).lon(minid);
                        front_lat = rhb_group.(DN).lat(minid);
                        seg(n).front_loc = [front_lon, front_lat];
                        seg(n).wavelength =  feature_wvlen;                        % take this from the scalogram;
                        
                        % setment 02: peak -> trough 2
                        xmask = (traj>=pks.loc(ip)) & (traj<=xlocs_sorted(2));
                        xmask2 = (dist_equal>=pks.loc(ip)) & (dist_equal<=xlocs_sorted(2));
                        
                        for iv = 1:length(VOIs)
                            varn = VOIs{iv};
                            seg(n+1).(varn) = rhb_group.(DN).(varn)(xmask);
                        end
                        seg(n+1).SSTgrad = sst_grad_abs(2) .* sign(mode(ship_wind_align(xmask2)));
                        seg(n+1).front_loc_traj = 0.5*(xlocs_sorted(2)+ pks.loc(ip));
                        % translate back to the lat-lon coordinate.
                        [minval, minid] = min(abs(traj - 0.5*(xlocs_sorted(2)+ pks.loc(ip))));
                        front_lon = rhb_group.(DN).lon(minid);
                        front_lat = rhb_group.(DN).lat(minid);
                        seg(n+1).front_loc = [front_lon, front_lat];
                        seg(n+1).wavelength = feature_wvlen;                      % take this from the scalogram;
                        
                        cnt = cnt +2;
                        
                        figure(10)
                        subplot(2,2,[2,4])
                        hold on;
                        if sign(seg(n).SSTgrad) >0      %
                            c1 = 'r';
                            c2 = 'b';
                        else
                            c1 = 'b';
                            c2='r';
                        end
                        plot(seg(n).lon, seg(n).lat,'o','color',c1);
                        hold on
                        plot(seg(n).front_loc(1), seg(n).front_loc(2),'*k','linewidth',1.2, 'markersize',10);
                        
                        plot(seg(n+1).lon, seg(n+1).lat,'o','color',c2);
                        plot(seg(n+1).front_loc(1), seg(n+1).front_loc(2),'*k','linewidth',1.2, 'markersize',10);
                        
                        text(seg(n).front_loc(1)+0.05, seg(n).front_loc(2), num2str(n));
                        text(seg(n+1).front_loc(1)+0.05, seg(n+1).front_loc(2), num2str(n+1));
                        
                        
                        %
                        figure(10);
                        subplot(2,2,1);
                        text(seg(n).front_loc_traj, pks.val(ip)-abs(sst_dif(1))/2 + 0.04, num2str(n),'fontsize',12, 'fontweight','bold');
                        text(seg(n+1).front_loc_traj, pks.val(ip)-abs(sst_dif(2))/2+0.04, num2str(n+1),'fontsize',12, 'fontweight','bold');
                        
                        
                        
                        
                        %% if gradient is only taken from one side:
                    elseif length(xlocs)==1
                        if xlocs<pks.loc(ip)
                            xmask = (traj>=xlocs) & (traj<=pks.loc(ip));
                            xmask2 = (dist_equal>=xlocs) & (dist_equal<=pks.loc(ip));
                        else
                            xmask = (traj>=pks.loc(ip)) & (traj<=xlocs);
                            xmask2 = (dist_equal>=pks.loc(ip)) & (dist_equal<=xlocs);
                        end
                        
                        for iv = 1:length(VOIs)
                            varn = VOIs{iv};
                            seg(n).(varn) = rhb_group.(DN).(varn)(xmask);
                        end
                        seg(n).SSTgrad = sst_grad_abs(1) .* sign(mode(ship_wind_align(xmask2)));
                        seg(n).wavelength = feature_wvlen;
                        
                        seg(n).front_loc_traj = 0.5*(xlocs+ pks.loc(ip));        % km;
                        % translate back to the lat-lon coordinate.
                        [minval, minid] = min(abs(traj - 0.5*(xlocs+ pks.loc(ip))));
                        front_lon = rhb_group.(DN).lon(minid);
                        front_lat = rhb_group.(DN).lat(minid);
                        seg(n).front_loc = [front_lon, front_lat];
                        
                        cnt = cnt +1;
                        
                        
                        figure(10);
                        subplot(2,2,[2,4])
                        hold on;
                        if sign(seg(n).SSTgrad) >0      %
                            c1 = 'r';
                            c2 = 'b';
                        else
                            c1 = 'b';
                            c2='r';
                        end
                        
                        hold on
                        plot(seg(n).lon, seg(n).lat,'o','color',c1);
                        plot(seg(n).front_loc(1), seg(n).front_loc(2),'*k','linewidth',1.2,'markersize',10);
                        text(seg(n).front_loc(1)+0.05, seg(n).front_loc(2), num2str(n));
                        
                        %
                        %                     if n==1
                        %                         plot(rhb_group.(DN).lon, rhb_group.(DN).lat,'--b');
                        %                         hold on
                        %                         % overlay the surface wind direction:
                        %                         wndir = get_cartesian_direction(rhb_group.(DN).wdir, 'Meteo');
                        %                         u = rhb_group.(DN).wspd_10N.* cosd(wndir);
                        %                         v = rhb_group.(DN).wspd_10N.* sind(wndir);
                        %
                        %                         %datetick(hb2,'x','mmmdd-hhZ','keepticks');
                        %                         %scatter(rhb_group.(DN).lon, rhb_group.(DN).lat, 10, rhb_group.(DN).local_time);
                        %                         scatter(rhb_group.(DN).lon, rhb_group.(DN).lat+0.25, 20, rhb_group.(DN).wspd_10N,'filled');
                        %                         hb2=colorbar;
                        %                         set(get(hb2,'xlabel'),'string','m/s');
                        %
                        %                         quiver(rhb_group.(DN).lon(1:6:end), rhb_group.(DN).lat(1:6:end)+0.25, ...
                        %                             u(1:6:end), v(1:6:end), 0.5,'k');
                        %                         axis('square');
                        %                     end
                        
                        
                        figure(10);
                        subplot(2,2,1);
                        text(seg(n).front_loc_traj, pks.val(ip)-sst_dif/2+0.04, num2str(n),'fontsize',12, 'fontweight','bold');
                        
                        
                        
                        
                        
                    end
                    
                    disp(['cumulative segment number:', num2str(cnt)]);
                    %pause
                    
                    
                end                                                                % end of manual peak selection
                
            end
            
            %% extra code dealing with a trough at the end of the records:
            
            if ip==num_pks
                
                % add an extra loop if the last trough is closer to the end of the SST
                % records:
                % find gradient between the last trough and the end point at feature
                % wavelength.
                sst_dif = [];
                sst_grad_abs =[];
                if (trghs.loc(end) > pks.loc(end))
                    % calcuate the feature wavelength associated with the last trough:
                    trwvlen = find_feature_wavelen(trghs.loc(end), dist_equal, mag_wt, WL);
                    
                    if (dist_equal(end)-trghs.loc(end))>trwvlen
                        % calculate the gradient at half feature wavelen:
                        xend = trghs.loc(end) + 0.75*trwvlen;
                        sst_end = interp1(dist_equal, SST_lowpassed, xend);
                        
                    else
                        % calculate gradient between the last trough and the end point
                        xend = dist_equal(end);
                        sst_end = SST_lowpassed(end);
                        
                    end
                    sst_dif(1) = (sst_end - (-trghs.val(end)));
                    sst_grad_abs(1) = sst_dif./(xend-trghs.loc(end));
                    
                    if flag_revisit_last_trgh
                        % compute the gradient to the left of the troughs
                        xlft = trghs.loc(end)-0.5*trwvlen;
                        sst_lft = interp1(dist_equal, SST_lowpassed, xlft);
                        sst_dif(2) = sst_lft-(-trghs.val(end));
                        sst_grad_abs(2) = sst_dif(2) ./ (xlft - trghs.loc(end));
                        xlocs = [xend, xlft];
                    else
                        xlocs = [xend];
                    end
                    
                    [xlocs_sorted, sxid] = sort(xlocs);
                    % sort the gradient and difference:
                    sst_dif = sst_dif(sxid);
                    sst_grad_abs = sst_grad_abs(sxid);
                    
                    figure(10);
                    subplot(2,2,1);
                    if length(xlocs)==2
                        plot(xlocs,[sst_end, sst_lft],'ob','linewidth',1.2,'markersize',13.5);
                    else
                        plot(xlocs,[sst_end],'ob','linewidth',1.2,'markersize',13.5);
                    end
                    
                    
                    n = cnt+1;
                    
                    ship_wind_align = rhb_group.(DN).ship_wind_align_mask_interp;
                    
                    
                    if length(xlocs_sorted)==2
                        % segment 01: xleft 1 -> trgh N
                        xmask = (traj>=xlocs_sorted(1)) & (traj<=trghs.loc(end));          % traj has the same length size as the time-based variables.
                        xmask2 = (dist_equal>=xlocs_sorted(1)) & (dist_equal<=trghs.loc(end));
                        for iv = 1:length(VOIs)
                            varn = VOIs{iv};
                            seg(n).(varn) = rhb_group.(DN).(varn)(xmask);
                        end
                        seg(n).SSTgrad = sst_grad_abs(1) .* mode(sign(ship_wind_align(xmask2)));
                        seg(n).front_loc_traj = 0.5*(xlocs_sorted(1)+ trghs.loc(end));        % km;
                        % translate back to the lat-lon coordinate.
                        [minval, minid] = min(abs(traj - 0.5*(xlocs_sorted(1)+ trghs.loc(end))));
                        front_lon = rhb_group.(DN).lon(minid);
                        front_lat = rhb_group.(DN).lat(minid);
                        seg(n).front_loc = [front_lon, front_lat];
                        seg(n).wavelength =  trwvlen;                        % take this from the scalogram;
                        
                        % setment 02: peak -> trough 2
                        xmask = (traj>=trghs.loc(end)) & (traj<=xlocs_sorted(2));
                        xmask2 = (dist_equal>=trghs.loc(end)) & (dist_equal<=xlocs_sorted(2));
                        for iv = 1:length(VOIs)
                            varn = VOIs{iv};
                            seg(n+1).(varn) = rhb_group.(DN).(varn)(xmask);
                        end
                        seg(n+1).SSTgrad = sst_grad_abs(2) .* sign(mode(ship_wind_align(xmask2)));
                        seg(n+1).front_loc_traj = 0.5*(xlocs_sorted(2)+ trghs.loc(end));
                        % translate back to the lat-lon coordinate.
                        [minval, minid] = min(abs(traj - seg(n+1).front_loc_traj));
                        front_lon = rhb_group.(DN).lon(minid);
                        front_lat = rhb_group.(DN).lat(minid);
                        seg(n+1).front_loc = [front_lon, front_lat];
                        seg(n+1).wavelength = trwvlen;                      % take this from the scalogram;
                        
                        cnt = cnt +2;
                        
                        figure(10)
                        subplot(2,2,[2,4])
                        hold on;
                        if sign(seg(n).SSTgrad) >0      %
                            c1 = 'r';
                            c2 = 'b';
                        else
                            c1 = 'b';
                            c2='r';
                        end
                        plot(seg(n).lon, seg(n).lat,'o','color',c1);
                        hold on
                        plot(seg(n).front_loc(1), seg(n).front_loc(2),'*k','linewidth',1.2, 'markersize',10);
                        
                        plot(seg(n+1).lon, seg(n+1).lat,'o','color',c2);
                        plot(seg(n+1).front_loc(1), seg(n+1).front_loc(2),'*k','linewidth',1.2, 'markersize',10);
                        
                        text(seg(n).front_loc(1)+0.05, seg(n).front_loc(2), num2str(n));
                        text(seg(n+1).front_loc(1)+0.05, seg(n+1).front_loc(2), num2str(n+1));
                        
                        
                        %
                        figure(10);
                        subplot(2,2,1);
                        text(seg(n).front_loc_traj, -trghs.val(end)+abs(sst_dif(1))/2 + 0.04, num2str(n),'fontsize',12, 'fontweight','bold');
                        text(seg(n+1).front_loc_traj, -trghs.val(end)+abs(sst_dif(2))/2+0.04, num2str(n+1),'fontsize',12, 'fontweight','bold');
                        
                        
                        
                        
                    elseif length(xlocs_sorted)==1
                        xmask = (traj>=trghs.loc(end)) & (traj<=xlocs_sorted);
                        xmask2 = (dist_equal>=trghs.loc(end)) & (dist_equal<=xlocs_sorted);
                        
                        for iv = 1:length(VOIs)
                            varn = VOIs{iv};
                            seg(n).(varn) = rhb_group.(DN).(varn)(xmask);
                        end
                        seg(n).SSTgrad = sst_grad_abs .* sign(mode(ship_wind_align(xmask2)));
                        seg(n).wavelength = trwvlen;
                        
                        seg(n).front_loc_traj = 0.5*(xend + trghs.loc(end));        % km;
                        % translate back to the lat-lon coordinate.
                        [minval, minid] = min(abs(traj -seg(n).front_loc_traj));
                        front_lon = rhb_group.(DN).lon(minid);
                        front_lat = rhb_group.(DN).lat(minid);
                        seg(n).front_loc = [front_lon, front_lat];
                        
                        % make plot on Figure 10 again:
                        figure(10);
                        subplot(2,2,1);
                        text(seg(n).front_loc_traj, -trghs.val(end)+abs(sst_dif)/2 + 0.04, num2str(n),'fontsize',12, 'fontweight','bold');
                        
                        
                        subplot(2,2,[2,4])
                        hold on;
                        if sign(seg(n).SSTgrad) >0      %
                            c1 = 'r';
                            c2 = 'b';
                        else
                            c1 = 'b';
                            c2='r';
                        end
                        
                        hold on
                        plot(seg(n).lon, seg(n).lat,'o','color',c1);
                        plot(seg(n).front_loc(1), seg(n).front_loc(2),'*k','linewidth',1.2,'markersize',10);
                        text(seg(n).front_loc(1)+0.05, seg(n).front_loc(2), num2str(n));
                    end
                    
                    disp(['cumulative segment number:', num2str(cnt)]);
                    
                end
                
                
                
                
                
                figure(10);
                subplot(2,2,[2,4]);
                hold on;
                plot(rhb_group.(DN).lon, rhb_group.(DN).lat,'--b');
                plot(rhb_group.(DN).lon(1), rhb_group.(DN).lat(1),'*k','markersize',13);
                plot(rhb_group.(DN).lon(end), rhb_group.(DN).lat(end),'dk','markersize',13);
                
                % overlay the surface wind direction:
                wndir = get_cartesian_direction(rhb_group.(DN).wdir, 'Meteo');
                u = rhb_group.(DN).wspd_10N.* cosd(wndir);
                v = rhb_group.(DN).wspd_10N.* sind(wndir);
                %scatter(rhb_group.(DN).lon, rhb_group.(DN).lat, 10, rhb_group.(DN).local_time);
                
                scatter(rhb_group.(DN).lon, rhb_group.(DN).lat+0.05, 20, rhb_group.(DN).wspd_10N,'filled');
                hb2=colorbar;
                set(get(hb2,'xlabel'),'string','m/s');
                %datetick(hb2,'x','mmmdd-hhZ','keepticks');
                caxis([median(rhb_group.(DN).wspd_10N)-2*std(rhb_group.(DN).wspd_10N), ...
                    median(rhb_group.(DN).wspd_10N)+2*std(rhb_group.(DN).wspd_10N)])
                quiver(rhb_group.(DN).lon(1:6:end), rhb_group.(DN).lat(1:6:end)+0.05, ...
                    u(1:6:end), v(1:6:end), 0.5,'k');
                axis('square')
                xlabel('longitude');
                ylabel('latitude');
                set(gca,'fontsize',14)
                
                
            end
            
        end                                                                    % end of peaks looping
        
    end
    
    
    
    if length(seg) > 0
        % save segments for each day:
        save([loc_datadir filesep DN filesep DN '_SST_front_segments.mat'],'seg');
    end
    figname = [DN '_SST_front_segments_selected_demo_improved_0913.jpg'];
    xc_savefig(hfig, [loc_datadir filesep DN], figname, [0 0 15 10]);
    
    pause
    
end                                                                        % end of the day loop.
 

   % Jan09: Amor
   % Jan18: Morse
   % Jan22: Morse: bad. need to improve the code's ability.
   % Jan23: Morse: similar issues as in Jan 22, the code doesn't find
   % distance to the left and right of the peak. (To-do: improve the code,
   % make a specific function to do this job) √
   % Jan24: Morse, interesting, wind perpenicular to the ship track, wind
   % do not feel the along-track SST differences as strong as it would if
   % the wind and ship are more aligned.
   
   % Feb04: similar issue as Jan23.
   % Feb08: worked well, seems that there are actually two eddies??  *
   % modify program a bit, so that if the first trough is considered. 
   
   % Feb11: There is issue with the end of record. 
   
   
%% now, let's find if there is a lead and lag in the wind and SST gradient. 
% try the standard lead-lag analysis; (cross correlation: xcorr)
% then also try this wavelet cross-correlation.
%
% r = xcorr(x,y,scaleopt);  % measure how much y lags or leads x;


% raw data:
wspd_10N =rhb_group.(DN).wspd_10N(moving_flag);
[c, lags] = xcorr(sst, wspd_10N,'normalized');
figure()
subplot(2,1,1);
yyaxis left
plot(dist_equal, U10N_interp,'-r','linewidth',1.2);
hold on
plot(traj, wspd_10N,'.r','linewidth',1.2);
plot(dist_equal, U10N_lowpassed, '--k','linewidth',1.2);
hold off

yyaxis right
plot(dist_equal, SST_lowpassed, '-');

subplot(2,1,2);
[cf, lags_f] = xcorr(SST_lowpassed, U10N_lowpassed,'normalized');

stem(lags_f, cf);



% try wavelet cross_correlation.
wt_sst = modwt(fliplr(SST_lowpassed));
wt_U10N = modwt(fliplr(U10N_interp));

[xc, xc_ci, lags] = modwtxcorr(wt_sst, wt_U10N);

lev = 1;
zerolag = floor(numel(xc{lev})/2+1);
%xlag = lags{lev}(zerolag-(170/2-1):zerolag+170/2).*2;
%xlag = lags{lev}.*2;
figure
for lev = 1:size(xc,1)
    subplot(size(xc,1),1,lev)
    zerolag = floor(numel(xc{lev})/2+1);
    %xlag = lags{lev}(zerolag-(140/2-1):zerolag+140/2).*2;
    xlag = lags{lev}.*2;
    
    %plot(xlag,xc{lev}(zerolag-(140/2-1):zerolag+140/2))
    plot(xlag, xc{lev});
    title(['Wavelet Cross-Correlation Sequence (level ' num2str(lev) ')'],['scale=' num2str(2^lev) 'to ' num2str(2^(lev+1)) 'km'])
    xlabel('lag (km)')
    ylabel('Cross-Correlation Coefficient')
    ylim([-0.8,0.8]);
    xlim([-500, 500]);
    
end
% 

% wavelet cross-correlation is used to measure similiarity between two
% signals at different scales. 
% wavelet cross-correlation is simply a scale-localized version of the
% usual cross-correlation between two signals. 

% an example:
t = 0:1/2000:1-1/2000;
x = sin(2*pi*200*t).*exp(-50*pi*(t-0.2).^2)+0.1*randn(size(t));
y = sin(2*pi*200*t).*exp(-50*pi*(t-0.5-(-0.3)).^2)+0.1*randn(size(t));

wx = modwt(x,'fk8',5);
wy = modwt(y,'fk8',5);

[xc,~,lags] = modwtxcorr(wx,wy,'fk8');
lev = 5;
zerolag = floor(numel(xc{lev})/2+1);
tlag = lags{lev}(zerolag-999:zerolag+1000).*(1/2000);
figure
plot(tlag,xc{lev}(zerolag-999:zerolag+1000))
title(['Wavelet Cross-Correlation Sequence (level '  num2str(lev) ')'])
xlabel('Time')
ylabel('Cross-Correlation Coefficient')

modwtcorr(wx,wy) ; % correlation coefficient at different scale. 