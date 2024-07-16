function [Tpick, Sxpick, Szpick, Gxpick, Gzpick, Offpick, topo, Errpick] = readpicks_pg(pickfile_GIMLI)

fidb = fopen(pickfile_GIMLI,'r'); % File opening %

npos = str2double(fgetl(fidb));
fgetl(fidb);

x = zeros(npos,1);
y = x; z = x;
for i = 1:npos
    test =  str2num(fgetl(fidb));
    x(i) = test(1);
    y(i) = test(2);
    z(i) = test(3);
end
if ~any(z) && any(y)
    z = y;
end

npicks = str2double(fgetl(fidb));
Tpick = zeros(npicks,1);
Sxpick = Tpick; Szpick = Tpick;
Gxpick = Tpick; Gzpick = Tpick;
Errpick = Tpick; Offpick = Tpick;

test = strsplit(fgetl(fidb));
ind_t = find(strcmp(test,'t')) - 1;
ind_s = find(strcmp(test,'s')) - 1;
ind_g = find(strcmp(test,'g')) - 1;
ind_err = find(strcmp(test,'err')) - 1;

for i = 1:npicks
    data = str2num(fgetl(fidb));
    Tpick(i) = data(ind_t) * 1000;
    Sxpick(i) = x(data(ind_s));
    Szpick(i) = z(data(ind_s));
    Gxpick(i) = x(data(ind_g));
    Gzpick(i) = z(data(ind_g));
    if ~isempty(ind_err)
        Errpick(i) = data(ind_err);
    end
    Offpick(i) = Gxpick(i) - Sxpick(i);
end
fclose(fidb);

ind_topo = find(diff(x) < 0,1);
if ~isempty(ind_topo)
    topo = [x(1:ind_topo) z(1:ind_topo)];
else
    topo = [x z];
end