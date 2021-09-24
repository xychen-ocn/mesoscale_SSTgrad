% wavelet analysis example:

load sumsin
figure()
plot(sumsin)
title('signal')

% perform a 3-level wavelet decomposition of the signal using the order 2
% Daubechies wavelet
[c,l] = wavedec(sumsin, 5 , 'db2');
% extract the coarse scale approx. coefficients and the detail coeff. from
% the decomposition.
approx = appcoef(c,l,'db2');
[cds] = detcoef(c,l,[1:5]);


% plot the coefficients:
figure()
subplot(4,1,1)
plot(approx)
title('Approximation Coefficients')

for i = 1:5
    subplot(5,1,i)
    plot(cds{i})
    title(['Level' num2str(i)])
end

subplot(4,1,3)
plot(cd2)
title('Level2')

subplot(4,1,4)
plot(cd1)
title('Level1')


test = rhb_group.Jan09.sst_seasnake;
dist_equal = [0:2:round(max(rhb_group.Jan09.traj))];
SST_interp = interp1(rhb_group.Jan09.traj, test ,dist_equal,'linear');

nmax = fix(log2(length(test)));
[c,l] = wavedec(test, 4, 'db2');
approx = appcoef(c,l,'db2');
cds = detcoef(c,l,[1:4]);

% plot the coefficients:
figure(1)
subplot(5,1,1)
plot(approx)
title('Approximation Coefficients')

for i = 1:4
    subplot(5,1,i+1)
    plot(cds{i})
    title(['Level' num2str(i)])
end

%
figure()
plot(test)
title('signal')

fs = 1/600;
[wt,f]=cwt(test, 'morse',fs);
[wt,fn]=cwt(detrend(detrend(SST_interp,1)),'amor');
% what is the magnitude of the sample?
cwt(detrend(detrend(SST_interp,1)),'amor');

trend_linear = SST_interp - detrend(SST_interp,1);

figure()
subplot(2,2,1)
yyaxis left
plot(SST_interp)
hold on
plot(trend_linear);
xlabel('# of samples');
ylabel('Celcius');

yyaxis right
plot(detrend(SST_interp,1));
hold on
ylabel('Celcius')


subplot(2,2,3)
% plot(dist_equal, SST_interp);
% hold on
% plot(dist_equal,SST_interp-trend_linear);
% hold on; 
plot(dist_equal, detrend(SST_interp,1));
title('detrended');
xlabel('km');
ylabel('Celcius');

% plot the scalogram: (magnitude of the wavelet transform)
subplot(2,2,[2,4])
%cwt(detrend(SST_interp),'amor');
wvlen = 1./fn .* 2;    % 1cycle/ (cycle per sample) * spatial resolution (km/sample) = km 
surface(dist_equal, wvlen , abs(wt));
shading flat
axis tight
xlabel('km');
ylabel('wavelengths (km)');
set(gca,'yscale','log');


%%
% 1. low-pass filter:
% wpass: normalized passband frequency: pi rad/sample.
% filter out wavelength smaller than 20km;
% wavenumber: 2*pi rad/20km;  resolution: 2km per sample
% rad per sample = wavenumber_cutoff * resolution 
wvlen_cutoff =20;   %units: km
wn_cutoff = 2/wvlen_cutoff ;  % pi rad/km
spatial_res = 2 ;   % 2km / sample
wpass = wn_cutoff * spatial_res;
SST_lowpassed = lowpass(SST_interp_dtr, wpass);  % retain frequency lower than the cutoff;


% use the location in x and length scale to compute SST gradient:
% perhaps I can also use the peak identification function in matlab in
% combination with the length scale. 
[pks, plocs, w, p] = findpeaks(SST_lowpassed, dist_equal, 'MinPeakDistance',20, 'MinPeakProminence',0.05);
[troughs, trlocs, w2, p2] = findpeaks(-SST_lowpassed, dist_equal, 'MinPeakDistance',20);


figure(10);clf;
subplot(2,1,1)
yyaxis left
plot(dist_equal,SST_interp)
hold on
plot(dist_equal,trend_linear);
xlabel('km');
ylabel('Celcius');
set(gca,'fontsize',14);


yyaxis right
plot(dist_equal,detrend(SST_interp,0),'o');
hold on
grid on
plot(plocs, pks,'*m','linewidth',1.2, 'markersize',12);
plot(trlocs, -troughs,'*b','linewidth',1.2, 'markersize',12);
plot(dist_equal, SST_lowpassed,'-','linewidth',1.2);
ylabel('Celcius')
set(gca,'fontsize',14);
xrange_p1 = get(gca,'xlim');
hold off
ylim([-0.6, 0.6]);


% plot the scalogram: (magnitude of the wavelet transform)
[wt,fn]=cwt(SST_lowpassed,'amor');
% use the location in space and wavelength to compute gradient for these
% features.
mag = abs(wt);
[DD, WL] = meshgrid(dist_equal, wvlen);
BW = imregionalmax(abs(wt));
% take the maximum that has a wavelength > 50km
C1 = WL>50;
C2 = mag>0.75*max(abs(wt(:)));
con = (C1&C2&BW);

% lengthscale of the local maximum :
wvlen_locmax = WL(con);
xloc_locmax = DD(con);

subplot(2,1,2)
%cwt(detrend(SST_interp),'amor');
levs = max(abs(wt(:))).*[0.5:0.1:1];
wvlen = 1./fn .* 2;    % 1cycle/ (cycle per sample) * spatial resolution (km/sample) = km 
pcolor(dist_equal, wvlen , abs(wt));
shading flat
hold on
plot(DD(con), WL(con),'om','MarkerSize',10);
[c,h] = contour(dist_equal, wvlen, abs(wt), levs,'k');
%clabel(c,h,levs);
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



%% find local maximum in the scalogram;
% use the scalogram to 

% now, estimate gradient:
% find pks that are cloest to the local maximum in the scalogram:
for i = 1:length(wvlen_locmax)
   [~, selid(i)]= min(abs(plocs-xloc_locmax(i)));
end
pk_sel.val = pks(selid);
pk_sel.loc = plocs(selid);
pk_sel.width = w(selid);
pk_sel.prominance = p(selid);


SST_interp_dtr= detrend(SST_interp,1);
dX = round(wvlen_locmax./2);  % half the wavelength
T_half_wvlen = zeros(length(selid),1);
dT = zeros(size(dX));
for ip = 1:length(selid)
    % find the SST at half wavelength from the peak:
    xs = [pk_sel.loc(ip)-dX(ip), pk_sel.loc(ip)+dX(ip)];
    tmp = interp1(dist_equal, SST_interp_dtr, xs);
    hold on
    plot(xs, tmp,'xk');
    dT_tmp = pk_sel.val(ip) - tmp;
    dT(ip) = max(abs(dT_tmp),[],'omitnan');
end
plot([xloc_locmax(1), xloc_locmax(1)+wvlen_locmax(1)],[-0.1, -0.1])
plot([xloc_locmax(2), xloc_locmax(2)-wvlen_locmax(2)],[-0.1, -0.1])
grad = dT./dX;


%% try polyfit (curve fitting the SST localized a given location) and then compute the gradient
% according to Li and Carbone (2012) paper:
section_mask = (dist_equal>=20)&(dist_equal<=66);
x0 =20;
xsub = dist_equal(section_mask)-40;
ysub = SST_lowpassed(section_mask);

error=[];
for i = 2:6
    [p{i}, S] = polyfit(xsub,ysub,i);
    yfit = polyval(p{i},xsub);
    error(i) = std(yfit - ysub);
    figure(12)
    plot(xsub, ysub,'.r');
    hold on
    plot(xsub, yfit,'-');
    hold off
    pause
end

%% instead of these methods, I will try calculating the gradient using the most standard method:
% compute the SST difference between consecutive peaks and troughs and then
% compute the gradient. (use manual selection to start with.)

% 


% The important data to save is 1: the peak --> the location of the maximum
% SST gradient sampled by the ship (--> set to half point between the
% maximum temperature and lowest temperature.)

% the time window, the distance traveled by the ship, 
