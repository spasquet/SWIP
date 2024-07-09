function [dpMIN,dpMAX]=param2mod(nameparam)

%%% S. Pasquet - V22.05.04
% Read .param dinver file
% [dpMIN,dpMAX]=param2mod(nameparam)

dir_param=fileparts(nameparam);
if strcmp(dir_param,'')==1
    dir_param=pwd;
end
nameparamnew=[nameparam(1:end-5),'tar'];

movefile(nameparam,nameparamnew,'f');
% [~,~]=unix_cmd(['mv -f ',nameparam_unix,' ',nameparamnew_unix]);

if ismac == 1
    unix_cmd(['gtar xf ',paramname,' contents.xml']);
elseif ispc == 1
    if ispc_wsl == 0
        [~,~]=unix_cmd(['tar -xf ',nameparamnew,' --force-local']);
        [~,~]=unix_cmd(['mv contents.xml ',dir_param]);
    else
        [~,~]=unix_cmd(['tar -xf ',nameparamnew,' -C ',dir_param]);
    end
else
    [~,~]=unix_cmd(['tar -xf ',nameparamnew,' -C ',dir_param]);
end
movefile(nameparamnew,nameparam,'f');
% [~,~]=unix_cmd(['mv -f ',nameparamnew_unix,' ',nameparam_unix]);

fid0=fopen(fullfile(dir_param,'contents.xml'),'r');
try
shortname=fgets(fid0); % Read line
catch
    keyboard
end
n=0;
while isempty(strfind(shortname,'<shortName>Vs</shortName>'))==1
    shortname=fgets(fid0); % Read line
end
for i=1:4
    fgets(fid0); % Read line
end
ligne=fgets(fid0); % Read line
while isempty(strfind(ligne,'<ParamLayer name="Vs'))~=1
    n=n+1;
    ligne=fgets(fid0); % Read line
    while isempty(strfind(ligne,'<dhMin'))==1
        ligne=fgets(fid0); % Read line
    end
    dhMin(n) = str2double(ligne(strfind(ligne,'<dhMin>')+7:strfind(ligne,'</dhMin>')-1));
    ligne=fgets(fid0); % Read line
    dhMax(n) = str2double(ligne(strfind(ligne,'<dhMax>')+7:strfind(ligne,'</dhMax>')-1));
    fgets(fid0); % Read line
    ligne=fgets(fid0); % Read line
end

dpMIN = sum(dhMin(1:end-1));
dpMAX = sum(dhMax(1:end-1));

fclose(fid0);

delete(fullfile(dir_param,'contents.xml'));

end
