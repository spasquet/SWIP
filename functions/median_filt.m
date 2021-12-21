function velout = median_filt(vel,npts,b1,b2,test)

%%% S. Pasquet - V16.11.18
% Median filter
% velout = median_filt(vel,npts,b1,b2)
if nargin == 2
    b1 = 1;
    b2 = length(vel);
end

if nargin < 5
    test = 0;
end

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
    if test == 0
        velout(i)=median(vel(i:i+npts-1));
    else
        % abs(median(vel(i:i+npts-1))-velout(i))>abs(median(vel(i:i+npts-1))-2*std(vel(i:i+npts-1)))
        if abs(median(vel(i:i+npts-1))-velout(i))>abs(0.75*std(vel(i:i+npts-1)))
            velout(i)=median(vel(i:i+npts-1));
        end
    end
end