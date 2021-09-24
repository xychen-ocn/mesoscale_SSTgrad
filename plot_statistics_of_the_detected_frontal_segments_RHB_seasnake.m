% Purpose: This script is used to show some statistics of the frontal
% segments found from the RHB (8 days) + (Jan 14, Jan 30, and Feb 11) add this 3 other days;
%

clear all; close all;

rhb_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
load([rhb_datadir filesep 'rhb_daily_grouped_10min_data_0909latest.mat']);


% load the segment file from RHB (8 strong SST variance days)
DOIs = rhbdates.strong_SSTvar;

spatial_scale = [];
SST_grad = [];
SST_anom = [];
seg_timewindow = [];
for i = 1:length(DOIs)
    DN = datestr(DOIs(i), 'mmmdd');
    segdir = [rhb_datadir filesep DN];
    segFN = [DN '_SST_front_segments.mat'];
    
    load([segdir filesep segFN]);
    
    % get the spatial scale, SST gradient of this segment:
    wvlen_all = unique([seg.wavelength]);        % this is the feature wavelength identify from the wavelet analysis with Morse base wavelet.
    spatial_scale = [spatial_scale, wvlen_all];
    
    SSTgrad_all = [seg.SSTgrad];
    SST_grad = [SST_grad, SSTgrad_all];
    
   % [seg.SST_anom]
    SST_anom_all = unique([seg.SST_anom]);           % (size of the SST anomly, the gradient is computed with half wavelength)
    SST_anom = [SST_anom, SST_anom_all];                    % there is potential to double count here..
    
    % --> get the survey time for the selected segment:
    seg_timewindow = [seg_timewindow, [seg.local_time]];
    
    radsonde_launch_timeUTC = DOIs(i)+2.75/24: 4/24 :DOIs(i)+1.25;
    
    hfig2= figure(3);
    set(hfig2, 'Name','radiosonde and front seg timing');
    subplot(5,2,i)
    % mark day and night zone:
    for iseg = 1:length(seg)
        yseg = 0.5*ones(size(seg(iseg).local_time));
        plot(seg(iseg).local_time, yseg, '-', 'linewidth',1.2);
        hold on;

        if std(seg(iseg).lon)>std(seg(iseg).lat)
            tx = interp1( seg(iseg).lon, seg(iseg).local_time, seg(iseg).front_loc(1));
        else
            tx = interp1( seg(iseg).lat, seg(iseg).local_time, seg(iseg).front_loc(2));
        end
        tx = mean(seg(iseg).local_time);
        plot(tx, 0.5,'*k', 'linewidth',0.9,'markersize',6);
    end
    yrad = repmat([0;1], 1, length(radsonde_launch_timeUTC));
    xrad = repmat(radsonde_launch_timeUTC-4/24, 2, 1);
    plot(xrad, yrad, '--','color','b');
    
    sunrise_loctime = DOIs(i)+7/24;
    sunset_loctime = DOIs(i)+17/24;
    mark_day_night_range(gca, sunrise_loctime, sunset_loctime);
    hold on

   hold off
    
    %set(gca,'xtick',[DOIs(i)+2:2/24:DOIs(i)+1.5]);
    ylim([0 1])
    xlim([DOIs(i), DOIs(i)+1]);
    set(gca,'xtick',[DOIs(i):2/24:DOIs(i)+1]);
    datetick('x','hh','keepticks');
    if i==9 || i==10
    xlabel('local time')
    end
    title(DN)

        
    % get mean, standard deviation on the neutral 10-m wind speed, fluxes
    % , MO_length over the segment.
    %for j = 1:length(seg)
    
    %end
    
    
    
    
    
end

seg_all.spatial_scale = spatial_scale;
seg_all.SST_grad = SST_grad;
seg_all.SST_anom = SST_anom;
seg_all.local_time = seg_timewindow;



% produce a histogram showing the range and distribution of the spatial scale of the SST gradients, 
% the magnitude of the SST gradients, the surface wind speed.
hfig1 = figure(1);
set(hfig1,'Name', 'histogram SST grad');

subplot(3,1,1)
histogram(abs(seg_all.SST_anom),[0:0.1:1],'Normalization','probability');
xlabel('SST anomaly (K)')
ylabel('Probability');
set(gca,'fontsize',14);


subplot(3,1,2)
histogram(seg_all.SST_grad,[-0.05:0.01:0.05], 'Normalization','probability');
xlabel('SST gradient (K/km)')
ylabel('Probability');
set(gca,'fontsize',14);


subplot(3,1,3)
histogram(seg_all.spatial_scale,[30:5:100],'Normalization','probability');
xlabel('Spatial scale of the SST gradient (km)');
ylabel('Probability');
set(gca,'fontsize',14);
ylim([0 0.3]);

figsvdir = [rhb_datadir filesep '../sst_variability_alongtrack/front_segments'];
xc_savefig(hfig1, figsvdir, ['statistics_of_selected_frontal_segments_on_10RHBdays_' num2str(length(seg_all.SST_grad)) 'segs.jpg'],[0,0, 10,8]);


% Collocate the frontal segments with radiosonde measurements to then obtain
% a histogram of the geostrophic wind speed (at 4km), and estimate of the
% ABL regimes (Zi/L, Zi is the inversion height, L is the Obukov length)

% we know that the radiosonde is lauched daily every 4 hours (6/day)
xc_savefig(hfig2, figsvdir, ['collocate_radiosonde_with_selected_fronts_time_map_10days.jpg'], [0, 0, 10, 12]);

