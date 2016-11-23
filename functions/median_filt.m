function velout = median_filt(vel,npts,b1,b2)

%%% S. Pasquet - V16.11.18
% Median filter
% velout = median_filt(vel,npts,b1,b2)

vel=vel';
velout=vel;
vel=[repmat(vel(1),fix(npts/2),1);vel;repmat(vel(end),fix(npts/2),1)];
indstart=find(isnan(vel)~=1,1,'first');
indend=find(isnan(vel)~=1,1,'last');
if indstart>1
    vel(indstart-fix(npts/2):indstart-1)=vel(indstart);
end
if indend<length(vel)
    vel(indend+1:indend+fix(npts/2))=vel(indend);
end
for i=b1:b2
    velout(i)=median(vel(i:i+npts-1));
end