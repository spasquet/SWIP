% function batch_rename(dir_dat)

%%% S. Pasquet - V17.09.05
% Batch rename files for seg2su or ascii2su processing
% batch_rename(dir_dat)

if exist('dir_dat','var')==0 || isempty(dir_dat)==1
    dir_dat=uigetdir('./','Select folder containing files to rename');
    if dir_dat==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    datstruct=[dir(fullfile(dir_dat,'*txt'));dir(fullfile(dir_dat,'*dat'));dir(fullfile(dir_dat,'*su'));dir(fullfile(dir_dat,'*sgy'))];
end

base_file_no=1000;

for ii=1:length(datstruct)
    [folder,nametarg,ext]=fileparts(datstruct(ii).name);
    for i=1:length(nametarg)
        test=str2double(nametarg(i));
        if isnan(test)==0 && isreal(test)==1
            file_no=str2double(nametarg(i:end));
            filename=fullfile(dir_dat,[num2str(base_file_no+file_no),ext]);
            movefile(fullfile(dir_dat,datstruct(ii).name),filename);
            break
        end
    end
end