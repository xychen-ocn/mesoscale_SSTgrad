function seg = find_fronts_and_compute_gradient(rhb_daily, moving_flag, dist_equal, SST_lowpassed, pks, trghs, mag_wt, WL,traj, VOIs, hfig)
% Purpose: detect SST fronts at selected scales by finding the slope around
%          a detected peak that has a wavelength above a threshold.
% return: the SST front segments and the associated gradients.
%
% Inputs: 
%
% Outputs:
% This function contains 2 nested fuctions.

%global pks traj dist_equal VOIs rhb_daily

cnt = 0;
num_pks = length(pks.val);
%min_dist = 10 ;     % units: km (hard coded) 12.5
seg = struct([]);                                                          % initialize the segment structure.
flag_revisit_trough = false;                                               % initialize trough.
thres_wvlen = 20;                                                         % 20km 
min_dist = 0.5*thres_wvlen; 

%moving_flag = rhb_daily.in_transit_flag;
rhb_lon = rhb_daily.lon(moving_flag);
rhb_lat = rhb_daily.lat(moving_flag);

for ip = 1:num_pks
    
    
    if ip == 1
        if (trghs.loc(1) < pks.loc(1)) & (abs(-trghs.val(1)-SST_lowpassed(1))>0.15)
            % if the beginning of the record is actually closer to
            % the first trough:
            dx = trghs.loc(1) - dist_equal(1);
            trwvlen = find_feature_wavelen(trghs.loc(1), dist_equal, mag_wt, WL);
            if trwvlen > thres_wvlen
                first_trgh.loc = trghs.loc(1);
                first_trgh.val = -trghs.val(1);
                first_trgh.trwvlen = trwvlen;
                frac = 0.5;
                if dx > frac*trwvlen
                    sidestr = 'left';
                    seg_tmp=compute_gradient_oneside_of_trough(first_trgh, dist_equal, SST_lowpassed,sidestr, frac);
                    seg = [seg, seg_tmp(end)];
                    disp(['cumulative segment number:', num2str(cnt)]);
                end
            end
        end
    end
    
    % get the mag_wt at the xlocs of the pks:
    xid = find(dist_equal==pks.loc(ip));
    TF = islocalmax(mag_wt(:,xid));
    [feature_wvlen, maxid] = max(WL(TF,xid));
    %[maxval, maxid] = max(mag_wt(:,xid));
    %feature_wvlen = WL(maxid,xid);
    %feature_wvlen=find_feature_wavelen(pks.loc(ip), dist_equal, mag_wt, WL);    
    good_peaks = (max(mag_wt(TF,xid))>0.70*max(mag_wt(:))) & (round(feature_wvlen)>=thres_wvlen);
    
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
            
            %% find fronts near a peak, look for it to both to the left and right of the peak:
            % ---  i. find the consecutive troughs, the cloestest one next to the peak from the left and right.
            dx =  trghs.loc - pks.loc(ip);
            
            crit = dx(1:end-1).*dx(2:end);
            iloc = find(crit<0);
            
            if ~isempty(iloc)
            tr_ids0 =[iloc, iloc+1];
            dx_sel0 = dx(iloc:iloc+1);                                                % initial selection
            
            
            % test if the selected troughs satisfy some criteria in distance:
            % go through each selected trough:
            trcnt =0;
            trs_xloc=[];
            trs_sst=[];
            for i = 1:2
                dx_upperlimit = abs(dx_sel0(i))<=feature_wvlen;
                dx_lowerlimit = abs(dx_sel0(i))>min_dist;
                cond_dist = dx_lowerlimit & dx_upperlimit;
                
                if cond_dist       % if this is true:
                    trcnt = trcnt+1;
                    trs_xloc(trcnt) = trghs.loc(tr_ids0(i));
                    trs_sst(trcnt) = -trghs.val(tr_ids0(i));
                    
                else
                    if ~dx_lowerlimit      % less than lower limit:
                        % two possible situations:
                        % a. the next peak is in between this trough and its
                        % consecutive trough.
                        
                        if ip<num_pks
                            
                            if dx_sel0(i)>0                                % trough to the right of the peak
                                nxtid = 1;
                            else                                           % trough to the left of the peak
                                nxtid = -1;
                            end
                            next_peak_loc = pks.loc(ip+nxtid);                 % or -1...
                            if length(trghs.loc)>tr_ids0(i)
                            consecutive_trgh_loc = trghs.loc(tr_ids0(i)+nxtid);
                            %if next_peak_loc < consecutive_trgh_loc       % this condition is not so good;; corrected
                            if (next_peak_loc - consecutive_trgh_loc)*(next_peak_loc-trghs.loc(tr_ids0(i))) < 0
                                % -->(next peak loc is within this trough and the conseutive trough
                                % do nothing, do not count this current trough as a valid
                                % option.
                                continue
                            else  % next_peak_loc > consecutive_trgh_loc:
                                % calculate gradient between the current peak and the
                                % consecutive trghs:
                                dx_tmp = abs(consecutive_trgh_loc - pks.loc(ip));
                                if dx_tmp < feature_wvlen
                                    trcnt = trcnt+1;
                                    trs_xloc(trcnt) = consecutive_trgh_loc;
                                    trs_sst(trcnt) =  -trghs.val(tr_ids0(i)+nxtid);    %+1 or -1
                                else
                                    % interpolation at half feature wavelength to find SST
                                    % gradient.
                                    trcnt = trcnt+1;
                                    trs_xloc(trcnt) = pks.loc(ip) + sign(nxtid).*0.5*feature_wvlen;
                                    trs_sst(trcnt) = interp1(dist_equal, SST_lowpassed, trs_xloc(trcnt),'linear');
                                    
                                    % need to revisit this consecutive troughs at the end
                                    % to compute the gradient at its half feature
                                    % wavelength.
                                    
                                    % find the feature wavelength for the trough:
                                    trwvlen = find_feature_wavelen(consecutive_trgh_loc, dist_equal, mag_wt, WL);
                                    if trwvlen > thres_wvlen
                                        dist_remained = dx_tmp - 0.5*feature_wvlen;
                                        if dist_remained>=0.5*trwvlen
                                            flag_revisit_trough = true;
                                            revisit_trgh.loc = consecutive_trgh_loc;
                                            revisit_trgh.trwvlen= trwvlen;
                                            revisit_trgh.val = -trghs.val(tr_ids0(i)+nxtid);  % or -1
                                            
                                            if nxtid==1
                                                revisit_trgh.side = 'left';
                                            else
                                                revisit_trgh.side = 'right';
                                            end
                                            
                                        end
                                    end
                                    
                                end
                            end
                            else
                                continue
                            end
                        end
                        
                        % if ip==num_pks, and less than lower limit, than the trough will
                        % be skipped -> do nothing.
                        
                        
                    elseif ~dx_upperlimit     %larger than feature wavelength:
                        % interpolation:
                        trcnt = trcnt+1;
                        if dx_sel0(i)>0
                            nxtid = 1;
                        else
                            nxtid = -1;
                        end
                        trs_xloc(trcnt) = pks.loc(ip) + sign(nxtid).* 0.5*feature_wvlen;       % this is the problem, could also be less than 0.5WL (on the other side.)
                        trs_sst(trcnt) = interp1(dist_equal, SST_lowpassed, trs_xloc(trcnt), 'linear');
                        
                        % find the feature wavelength for the trough:
                        trwvlen = find_feature_wavelen(trghs.loc(tr_ids0(i)), dist_equal, mag_wt, WL);
                        if trwvlen > thres_wvlen
                            dist_remained = abs(dx_sel0(i)) - 0.5*feature_wvlen;
                            if dist_remained>0.5*trwvlen
                                flag_revisit_trough = true;
                                revisit_trgh.loc = trghs.loc(tr_ids0(i));
                                revisit_trgh.trwvlen= trwvlen;
                                revisit_trgh.val = -trghs.val(tr_ids0(i));
                                if nxtid==1
                                    revisit_trgh.side = 'left';
                                else
                                    revisit_trgh.side = 'right';
                                end
                            end
                        end
                        
                        
                        
                    end
                end
            end
            
            %  show the frontal section in a figure. (sanity check);
%             figure(10)
%             subplot(2,2,1)
%             hold on
%             plot(trs_xloc, trs_sst,'ob','linewidth',1.5, 'markersize',13.5);
            
            % note: trs_xloc is already in ascending order.
            
            %% use the found SST fronts (conditional--> satisfy certain criteria) to compute SST gradient.
            % the code packet in the following section can be written in
            % the a function:
            
            seg_tmp = compute_gradient_bothside_of_peak(trs_xloc, trs_sst, feature_wvlen);
            if length(trs_xloc)==2
                seg = [seg, seg_tmp(end-1:end)];
            else
                seg=[seg,seg_tmp(end)];
            end
            
            %% code here deals with the troughs that needs to be revisited:
            if flag_revisit_trough
                
                % compute the gradient for one side:
                sidestr = revisit_trgh.side;
                seg_tmp=compute_gradient_oneside_of_trough(revisit_trgh, dist_equal, SST_lowpassed,sidestr);
                seg = [seg, seg_tmp(end)];
                disp(['cumulative segment number:', num2str(cnt)]);

            end

            
            else % peak is not within two troughs: This situation occur either at the beginning or the end of the SST records:
                
                disp(['ip=' num2str(ip), '; num_pks=', num2str(num_pks)]);
                % determin if the first SST record to the peak is > 0.5
                % wavelen, if so, two side, same if the peak is close to
                % the end of the record.
                
                % find the appropriate trough on a single side:
                dx_abs_ascending = sort(abs(dx));
                
                % take the first two (in case the minimal distance is
                % too short)
                dx_tmp = dx_abs_ascending(1:2);
                
                % see which one satisfy the condition:
                dx_upperlimit = dx_abs_ascending(1:2)<=feature_wvlen;
                dx_lowerlimit = dx_abs_ascending(1:2)>min_dist;
                cond_dist = dx_lowerlimit & dx_upperlimit;
                
                min_trgh_dist = min(dx_tmp(cond_dist));
                
                if ~isempty(min_trgh_dist)
                    % find location of this trough:
                    trs_xloc=[]; trs_sst=[];
                    iloc = find(abs(dx) == min_trgh_dist);
                    trs_xloc(1) = trghs.loc(iloc);
                    trs_sst(1) = -trghs.val(iloc);
                    
                else
                    continue
                end
                
                
                if ip==1
                    
                    % determine if the SST records to the left of the peak
                    % is long enough to compute gradient:
                   % if pks.loc(1) > trghs.loc(1)
                        dx = pks.loc(ip) - dist_equal(1);
                        if dx > 0.5*feature_wvlen
                            trs_xloc(2) = pks.loc(ip) - 0.5*feature_wvlen;
                            trs_sst(2) = interp1(dist_equal, SST_lowpassed, trs_xloc(2),'linear');
                        end
                        
                        % sort:
                        if length(trs_xloc)==2
                            [trs_xloc, sid] = sort(trs_xloc);
                            trs_sst = trs_sst(sid);
                        end
                        
                        seg_tmp=compute_gradient_bothside_of_peak(trs_xloc, trs_sst, feature_wvlen);
                        if length(trs_xloc)==2
                            seg = [seg, seg_tmp(end-1:end)];
                        else
                            seg=[seg, seg_tmp(end)];
                        end
                        
                        
                   
                    
                    % need to deal with peak that is close to the beginning and
                    % end of the SST records:
                elseif ip==num_pks
                    
                    % determine if the SST records to the right of the peak
                    % is long enough to compute gradient:                    

                    dx = abs(pks.loc(ip) - dist_equal(end));
                    if dx > 0.5*feature_wvlen
                        trs_xloc(2) = pks.loc(ip) + 0.5*feature_wvlen;
                        trs_sst(2) = interp1(dist_equal, SST_lowpassed, trs_xloc(2),'linear');
                    end
                    
                    % sort:
                    if length(trs_xloc)==2
                        [trs_xloc, sid] = sort(trs_xloc);
                        trs_sst = trs_sst(sid);
                    end
                    
                    seg_tmp=compute_gradient_bothside_of_peak(trs_xloc, trs_sst, feature_wvlen);
                    if length(trs_xloc)==2
                        seg = [seg, seg_tmp(end-1:end)];
                    else
                        seg=[seg, seg_tmp(end)];
                    end
                    
           
                    %
                end
                
                
                
            end
        end
    
  
    
        
    
    if ip == num_pks    % at the last peak:
        %% code that deal with records with the last trough that is closer to the end of the SST records:
        if (trghs.loc(end) > pks.loc(end))
            disp('last trough closer to the end');

            trwvlen = find_feature_wavelen(trghs.loc(end), dist_equal, mag_wt, WL);
            disp(['trwvlen:', num2str(trwvlen)]);
            frac = 0.5;
            if trwvlen > thres_wvlen
                if (dist_equal(end)-trghs.loc(end))>frac*trwvlen
                    % calculate the gradient at half feature wavelen:
                    last_trgh.loc = trghs.loc(end);
                    last_trgh.trwvlen = trwvlen;
                    last_trgh.val = -trghs.val(end);
                    
                    seg_tmp=compute_gradient_oneside_of_trough(last_trgh, dist_equal, SST_lowpassed,'right',frac);
                    seg = [seg, seg_tmp(end)];
                    disp(['cumulative segment number:', num2str(cnt)]);
                end
            end
        end

        
    end
    
    
 
    end                                                                    % close condition: if it is a good peak
    
end                                                                        % end of peaks looping

disp(['--> total number of segments:' num2str(length(seg))]);
% %seg_tmp.SSTgrad
% seg.SSTgrad




%% nested functions:
function  seg=compute_gradient_oneside_of_trough(trghinfo,dist, SST_lowpassed, sidestr, frac)
switch nargin
    case 4
        wvlenfrac = 0.5;
    case 5
        wvlenfrac = frac;
end
n = cnt+1;
% find the SST at 0.5*trwvlen left from the trough of intersts
if strcmp(sidestr,'left')
    xloc = trghinfo.loc - wvlenfrac*trghinfo.trwvlen;
else
    xloc = trghinfo.loc + wvlenfrac*trghinfo.trwvlen;
end
sst = interp1(dist, SST_lowpassed, xloc);

if strcmp(sidestr,'left')
    sst_grad0 = (trghinfo.val - sst)./(wvlenfrac*trghinfo.trwvlen);
else
    sst_grad0 = (trghinfo.val - sst)./(-wvlenfrac*trghinfo.trwvlen);
end

% figure out how the ship is travel relative to the wind
% from trough to peak and from peak to troughs.
ship_wind_align = rhb_daily.ship_wind_align_mask_interp;


% relative to the ship:
if strcmp(sidestr,'left')
    xmask = (traj>=xloc) & (traj<=trghinfo.loc);
    xmask2 = (dist>=xloc) & (dist_equal<=trghinfo.loc);
else
    xmask = (traj>=trghinfo.loc) & (traj<=xloc);
    xmask2 = (dist>=trghinfo.loc) & (dist_equal<=xloc);
    
end

for iv = 1:length(VOIs)
    varn = VOIs{iv};
    moving_tmp  = rhb_daily.(varn)(moving_flag);
    seg(n).(varn) = moving_tmp(xmask);
end
seg(n).SSTgrad = sst_grad0 .* sign(mode(ship_wind_align(xmask2)));
seg(n).wavelength = trghinfo.trwvlen;

seg(n).front_loc_traj = 0.5*(xloc + trghinfo.loc);        % km;
% translate back to the lat-lon coordinate.
[minval, minid] = min(abs(traj -seg(n).front_loc_traj));
front_lon = rhb_lon(minid);
front_lat = rhb_lat(minid);
seg(n).front_loc = [front_lon, front_lat];
seg(n).SST_anom = (trghinfo.val - sst);

% make plot on Figure 10 again:
figure(10);
subplot(2,2,1);
hold on
plot(xloc, sst,'ob','linewidth',1.5, 'markersize',13.5);
text(seg(n).front_loc_traj, trghinfo.val+abs(trghinfo.val - sst)/2 + 0.04, num2str(n),'fontsize',12, 'fontweight','bold');


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
text(seg(n).front_loc(1), seg(n).front_loc(2)+0.025, num2str(n));


cnt = cnt +1;
end


function seg = compute_gradient_bothside_of_peak(trs_xloc, trs_sst, feature_wvlen)
%  show the frontal section in a figure. (sanity check);
figure(10)
subplot(2,2,1)
hold on
plot(trs_xloc, trs_sst,'ob','linewidth',1.5, 'markersize',13.5);

% note: trs_xloc is already in ascending order.

%% use the found SST fronts (conditional--> satisfy certain criteria) to compute SST gradient.
% now compute the SST gradient:
sst_dif = pks.val(ip) - trs_sst;
sst_grad_rhb = sst_dif ./ (pks.loc(ip) - trs_xloc) ;                      % this is sst_grad has sign associated with it, and it relative to the rhb ship direction.


% figure out how the ship is travel relative to the wind
% from trough to peak and from peak to troughs.
ship_wind_align = rhb_daily.ship_wind_align_mask_interp;


%% if gradient is taken from two sides of the peak
if length(trs_xloc)==2
    % segment 01: trough 1 -> peak0
    xmask = (traj>=trs_xloc(1)) & (traj<=pks.loc(ip));          % traj has the same length size as the time-based variables.
    for iv = 1:length(VOIs)
        varn = VOIs{iv};
        moving_tmp  = rhb_daily.(varn)(moving_flag);
        seg(n).(varn) = moving_tmp(xmask);
        %seg(n).(varn) = rhb_daily.(varn)(xmask);
    end
    % ===== a different mask that has the same length as dist_equal) === %
    xmask2 = (dist_equal>=trs_xloc(1)) & (dist_equal<=pks.loc(ip));
    seg(n).SSTgrad = sst_grad_rhb(1) .* mode(sign(ship_wind_align(xmask2)));
    seg(n).front_loc_traj = 0.5*(trs_xloc(1)+ pks.loc(ip));        % km;
    % translate back to the lat-lon coordinate.
    [minval, minid] = min(abs(traj - seg(n).front_loc_traj ));
    front_lon = rhb_lon(minid);
    front_lat = rhb_lat(minid);
    seg(n).front_loc = [front_lon, front_lat];
    seg(n).wavelength =  feature_wvlen;                        % take this from the scalogram;
    seg(n).SST_anom = sst_dif(1);
    
    % setment 02: peak -> trough 2
    xmask = (traj>=pks.loc(ip)) & (traj<=trs_xloc(2));
    for iv = 1:length(VOIs)
        varn = VOIs{iv};
        moving_tmp  = rhb_daily.(varn)(moving_flag);
        seg(n+1).(varn) = moving_tmp(xmask);
        %seg(n+1).(varn) = rhb_daily.(varn)(xmask);
    end
    
    xmask2 = (dist_equal>=pks.loc(ip)) & (dist_equal<=trs_xloc(2));
    seg(n+1).SSTgrad = sst_grad_rhb(2) .* sign(mode(ship_wind_align(xmask2)));
    
    seg(n+1).front_loc_traj = 0.5*(trs_xloc(2)+ pks.loc(ip));
    % translate back to the lat-lon coordinate.
    [minval, minid] = min(abs(traj - seg(n+1).front_loc_traj));
    front_lon = rhb_lon(minid);
    front_lat = rhb_lat(minid);
    seg(n+1).front_loc = [front_lon, front_lat];
    seg(n+1).wavelength = feature_wvlen;                      % take this from the scalogram;
    seg(n+1).SST_anom = sst_dif(2);

    
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
    
    text(seg(n).front_loc(1), seg(n).front_loc(2)+0.025, num2str(n));
    text(seg(n+1).front_loc(1), seg(n+1).front_loc(2)+0.025, num2str(n+1));
    
    
    %
    figure(10);
    subplot(2,2,1);
    text(seg(n).front_loc_traj, pks.val(ip)-abs(sst_dif(1))/2 + 0.04, num2str(n),'fontsize',12, 'fontweight','bold');
    text(seg(n+1).front_loc_traj, pks.val(ip)-abs(sst_dif(2))/2+0.04, num2str(n+1),'fontsize',12, 'fontweight','bold');
    
    
    
    
    %% if gradient is only taken from one side:
elseif length(trs_xloc)==1
    if trs_xloc<pks.loc(ip)
        xmask = (traj>=trs_xloc) & (traj<=pks.loc(ip));
        xmask2 = (dist_equal>=trs_xloc) & (dist_equal<=pks.loc(ip));
    else
        xmask = (traj>=pks.loc(ip)) & (traj<=trs_xloc);
        xmask2 = (dist_equal>=pks.loc(ip)) & (dist_equal<=trs_xloc);
    end
    
    for iv = 1:length(VOIs)
        varn = VOIs{iv};
        moving_tmp  = rhb_daily.(varn)(moving_flag);
        seg(n).(varn) = moving_tmp(xmask);
        %seg(n).(varn) = rhb_daily.(varn)(xmask);
    end
    seg(n).SSTgrad = sst_grad_rhb(1) .* sign(mode(ship_wind_align(xmask2)));
    seg(n).wavelength = feature_wvlen;
    
    seg(n).front_loc_traj = 0.5*(trs_xloc+ pks.loc(ip));        % km;
    % translate back to the lat-lon coordinate.
    [minval, minid] = min(abs(traj - seg(n).front_loc_traj ));
    front_lon = rhb_lon(minid);
    front_lat = rhb_lat(minid);
    seg(n).front_loc = [front_lon, front_lat];
    seg(n).SST_anom = sst_dif;                                             % knowing this and the gradient, I can return the dx where the gradient is computed, if needed.
    
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
    text(seg(n).front_loc(1), seg(n).front_loc(2)+0.025, num2str(n));
    
    figure(10);
    subplot(2,2,1);
    text(seg(n).front_loc_traj, pks.val(ip)-sst_dif/2+0.04, num2str(n),'fontsize',12, 'fontweight','bold');
    
    
    
    
    
end

disp(['cumulative segment number:', num2str(cnt)]);

end


end



%%%% obsolete:

% now compute the SST gradient:
%             sst_dif = pks.val(ip) - trs_sst;
%             sst_grad_rhb = sst_dif ./ (pks.loc(ip) - trs_xloc) ;                      % this is sst_grad has sign associated with it, and it relative to the rhb ship direction.
%             
%             
%             % figure out how the ship is travel relative to the wind
%             % from trough to peak and from peak to troughs.
%             ship_wind_align = rhb_daily.ship_wind_align_mask_interp;
%             
%             
%             %% if gradient is taken from two sides of the peak
%             if length(trs_xloc)==2
%                 % segment 01: trough 1 -> peak0
%                 xmask = (traj>=trs_xloc(1)) & (traj<=pks.loc(ip));          % traj has the same length size as the time-based variables.
%                 for iv = 1:length(VOIs)
%                     varn = VOIs{iv};
%                     seg(n).(varn) = rhb_daily.(varn)(xmask);
%                 end
%                 % ===== a different mask that has the same length as dist_equal) === %
%                 xmask2 = (dist_equal>=trs_xloc(1)) & (dist_equal<=pks.loc(ip));
%                 seg(n).SSTgrad = sst_grad_rhb(1) .* mode(sign(ship_wind_align(xmask2)));
%                 seg(n).front_loc_traj = 0.5*(trs_xloc(1)+ pks.loc(ip));        % km;
%                 % translate back to the lat-lon coordinate.
%                 [minval, minid] = min(abs(traj - seg(n).front_loc_traj ));
%                 front_lon = rhb_daily.lon(minid);
%                 front_lat = rhb_daily.lat(minid);
%                 seg(n).front_loc = [front_lon, front_lat];
%                 seg(n).wavelength =  feature_wvlen;                        % take this from the scalogram;
%                 
%                 % setment 02: peak -> trough 2
%                 xmask = (traj>=pks.loc(ip)) & (traj<=trs_xloc(2));
%                 for iv = 1:length(VOIs)
%                     varn = VOIs{iv};
%                     seg(n+1).(varn) = rhb_daily.(varn)(xmask);
%                 end
%                 
%                 xmask2 = (dist_equal>=pks.loc(ip)) & (dist_equal<=trs_xloc(2));
%                 seg(n+1).SSTgrad = sst_grad_rhb(2) .* sign(mode(ship_wind_align(xmask2)));
%                 
%                 seg(n+1).front_loc_traj = 0.5*(trs_xloc(2)+ pks.loc(ip));
%                 % translate back to the lat-lon coordinate.
%                 [minval, minid] = min(abs(traj - seg(n+1).front_loc_traj));
%                 front_lon = rhb_daily.lon(minid);
%                 front_lat = rhb_daily.lat(minid);
%                 seg(n+1).front_loc = [front_lon, front_lat];
%                 seg(n+1).wavelength = feature_wvlen;                      % take this from the scalogram;
%                 
%                 cnt = cnt +2;
%                 
%                 figure(10)
%                 subplot(2,2,[2,4])
%                 hold on;
%                 if sign(seg(n).SSTgrad) >0      %
%                     c1 = 'r';
%                     c2 = 'b';
%                 else
%                     c1 = 'b';
%                     c2='r';
%                 end
%                 plot(seg(n).lon, seg(n).lat,'o','color',c1);
%                 hold on
%                 plot(seg(n).front_loc(1), seg(n).front_loc(2),'*k','linewidth',1.2, 'markersize',10);
%                 
%                 plot(seg(n+1).lon, seg(n+1).lat,'o','color',c2);
%                 plot(seg(n+1).front_loc(1), seg(n+1).front_loc(2),'*k','linewidth',1.2, 'markersize',10);
%                 
%                 text(seg(n).front_loc(1)+0.05, seg(n).front_loc(2), num2str(n));
%                 text(seg(n+1).front_loc(1)+0.05, seg(n+1).front_loc(2), num2str(n+1));
%                 
%                 
%                 %
%                 figure(10);
%                 subplot(2,2,1);
%                 text(seg(n).front_loc_traj, pks.val(ip)-abs(sst_dif(1))/2 + 0.04, num2str(n),'fontsize',12, 'fontweight','bold');
%                 text(seg(n+1).front_loc_traj, pks.val(ip)-abs(sst_dif(2))/2+0.04, num2str(n+1),'fontsize',12, 'fontweight','bold');
%                 
%                 
%                 
%                 
%                 %% if gradient is only taken from one side:
%             elseif length(trs_xloc)==1
%                 if trs_xloc<pks.loc(ip)
%                     xmask = (traj>=trs_xloc) & (traj<=pks.loc(ip));
%                     xmask2 = (dist_equal>=trs_xloc) & (dist_equal<=pks.loc(ip));
%                 else
%                     xmask = (traj>=pks.loc(ip)) & (traj<=trs_xloc);
%                     xmask2 = (dist_equal>=pks.loc(ip)) & (dist_equal<=trs_xloc);
%                 end
%                 
%                 for iv = 1:length(VOIs)
%                     varn = VOIs{iv};
%                     seg(n).(varn) = rhb_daily.(varn)(xmask);
%                 end
%                 seg(n).SSTgrad = sst_grad_rhb(1) .* sign(mode(ship_wind_align(xmask2)));
%                 seg(n).wavelength = feature_wvlen;
%                 
%                 seg(n).front_loc_traj = 0.5*(xlocs+ pks.loc(ip));        % km;
%                 % translate back to the lat-lon coordinate.
%                 [minval, minid] = min(abs(traj - seg(n).front_loc_traj ));
%                 front_lon = rhb_daily.lon(minid);
%                 front_lat = rhb_daily.lat(minid);
%                 seg(n).front_loc = [front_lon, front_lat];
%                 
%                 cnt = cnt +1;
%                 
%                 
%                 figure(10);
%                 subplot(2,2,[2,4])
%                 hold on;
%                 if sign(seg(n).SSTgrad) >0      %
%                     c1 = 'r';
%                     c2 = 'b';
%                 else
%                     c1 = 'b';
%                     c2='r';
%                 end
%                 
%                 hold on
%                 plot(seg(n).lon, seg(n).lat,'o','color',c1);
%                 plot(seg(n).front_loc(1), seg(n).front_loc(2),'*k','linewidth',1.2,'markersize',10);
%                 text(seg(n).front_loc(1)+0.05, seg(n).front_loc(2), num2str(n));
%                 
%                 figure(10);
%                 subplot(2,2,1);
%                 text(seg(n).front_loc_traj, pks.val(ip)-sst_dif/2+0.04, num2str(n),'fontsize',12, 'fontweight','bold');
%                 
%                 
%                 
%                 
%                 
%             end
%             
%             disp(['cumulative segment number:', num2str(cnt)]);
            
