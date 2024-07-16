function [Tpick, Sxpick, Szpick, Gxpick, Gzpick, Offpick] = readpicks(input)

Tpick=[]; Sxpick=[]; Szpick=[]; Offpick=[]; Gxpick=[]; Gzpick=[];
Sx_sing = input(input(:,1) == 0,2);
nsrc = length(Sx_sing);

for i=1:nsrc;
    ind = find(input(:,1) == i);
    nrec = length(ind);
    I = input(:,1) == i;
    J = find(input(:,4) == i);
    Tpick=[Tpick;input(I,4)*1000];
    Sxpick=[Sxpick;repmat(input(J,2),nrec,1)];
    Szpick=[Szpick;repmat(input(J,3),nrec,1)];
    Gxpick=[Gxpick;input(I,2)];
    Gzpick=[Gzpick;input(I,3)];
    Offpick=[Offpick;input(I,2)-repmat(input(J,2),nrec,1)];
end