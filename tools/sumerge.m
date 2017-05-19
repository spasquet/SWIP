function sumerge(dir_su,filename_new,del_file)

%%% S. Pasquet - V17.05.10
% Merge multiple SU files
% sumerge(dir_su,filename_new,del_file)

if exist('dir_su','var')==0 || isempty(dir_su)==1
    dir_su=uigetdir('./','Select folder containing SU files');
    if dir_su==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a folder');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
    sustruct=dir(fullfile(dir_su,'*.su'));
end

if exist('filename_new','var')==0 || isempty(filename_new)==1
    filename_new = 'swip_profile_merged.su';
end

if exist('del_file','var')==0 || isempty(del_file)==1
    del_file = 0;
end

com1='cat';
for i=1:length(sustruct)
    com1 = [com1 ' ' sustruct(i).name];
end
com1=[com1 ' > ' filename_new];
unix(com1);

if del_file == 1
    delete(sustruct.name);
end