function feature_wvlen = find_feature_wavelen(xloc, dist_axis, mag_wt, WL)
xid = find(dist_axis ==xloc);
%TF = islocalmax(mag_wt(:,xid));
[maxval, maxid] = max(mag_wt(:,xid));
feature_wvlen = WL(maxid,xid);

return