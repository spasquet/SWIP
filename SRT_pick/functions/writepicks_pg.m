function writepicks_pg(pickfile_GIMLI,Tpick, Sxpick, Szpick, Gxpick, Gzpick, topo, Errpick)

% Save pick file for pyGIMLI
fidb=fopen(pickfile_GIMLI,'wt'); % File opening %
all_points = topo;

npoints = size(all_points,1);
fprintf(fidb,'%d\n',npoints);
fprintf(fidb,'%s\n','# x y z');
for ii = 1:npoints
    fprintf(fidb,'%f\t%f\t%f\n',all_points(ii,1),all_points(ii,2),0);
end
fprintf(fidb,'%d\n',length(Gxpick));
fprintf(fidb,'%s\n','# g s err t');
for jj = 1:length(Gxpick)
    [~,ind_s] = ismember([Sxpick(jj),Szpick(jj)],all_points,'rows');
    [~,ind_g] = ismember([Gxpick(jj),Gzpick(jj)],all_points,'rows');
    fprintf(fidb,'%d\t%d\t%f\t%f\n',ind_g,ind_s,Errpick(jj)/1000,(Tpick(jj))/1000);
end
fclose(fidb);
