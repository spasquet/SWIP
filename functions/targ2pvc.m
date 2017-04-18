function [freq,vel,deltac,modes]=targ2pvc(nametarg)

%%% S. Pasquet - V17.04.14
% Read .target dinver file
% [freq,vel,deltac,modes]=targ2pvc(nametarg)

dir_inv=fileparts(nametarg);
if strcmp(dir_inv,'')==1
    dir_inv=pwd;
end
copyfile(nametarg,nametarg(1:end-3));
if isunix==1
    [~,~]=unix(['tar -xf ',nametarg(1:end-3),' -C ',dir_inv]);
else
    [~,~]=unix(['tar -xf ',nametarg(1:end-3),' --force-local']);
    [~,~]=unix(['mv contents.xml ',dir_inv]);
end
delete(nametarg(1:end-3));

fid0=fopen(fullfile(dir_inv,'contents.xml'),'r');
for i=1:9
    fgets(fid0); % Read line
end
modalcurve=fgets(fid0); % Read line
n=0;
while strcmp(modalcurve(1:end-2),'<ModalCurve>')==1
    n=n+1;
    for i=1:2
        ligne=fgets(fid0); % Read line
    end
    nsamples=str2double(ligne(strfind(ligne,'<log>')+5:strfind(ligne,' samples')-1));
    for i=1:5
        ligne=fgets(fid0); % Read line
    end
    modes(n)=str2double(ligne(strfind(ligne,'<index>')+7:strfind(ligne,'</index>')-1));
    fgets(fid0); % Read line
    for i=1:nsamples
        fgets(fid0); % Read line
        ligne=fgets(fid0); % Read line
        freq{modes(n)+1}(i)=str2double(ligne(strfind(ligne,'<x>')+3:strfind(ligne,'</x>')-1));
        ligne=fgets(fid0); % Read line
        vel{modes(n)+1}(i)=1/str2double(ligne(strfind(ligne,'<mean>')+6:strfind(ligne,'</mean>')-1));
        ligne=fgets(fid0); % Read line
        deltaslow=str2double(ligne(strfind(ligne,'<stddev>')+8:strfind(ligne,'</stddev>')-1));
        deltac{modes(n)+1}(i)=deltaslow/((1/vel{modes(n)+1}(i))^2-deltaslow^2);
        for j=1:3
            fgets(fid0); % Read line
        end
    end
%     errorbar(freq{modes(n)+1},(vel{modes(n)+1}),deltac{modes(n)+1})
    fgets(fid0); % Read line
    modalcurve=fgets(fid0); % Read line
end
fclose(fid0);

delete(fullfile(dir_inv,'contents.xml'));

end
