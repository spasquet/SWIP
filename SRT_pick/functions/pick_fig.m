function [x, y, key] = pick_fig

%%% S. Pasquet - V17.04.20
% Get X,Y coordinates with mouse pick

k = waitforbuttonpress;
x=[]; y=[]; key='';
if k==0
    CP = get(gca, 'CurrentPoint');
    CP = CP(1,1:2);
    x  = CP(1);
    y  = CP(2);
elseif k==1
    key = get(gcf,'CurrentKey');
end
end