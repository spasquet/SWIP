function cm_saturation(sat)
%Changes saturation level of current figure use -1:1
ax=gca;
map=colormap(ax);

newmap=brighten(map,sat);


colormap(ax,newmap);