function [xnew,ynew] = convert_xzv(filevel,utm_start,utm_end)

%%% S. Pasquet - V17.05.19
% Convert 2D models into UTM coordinate system from X,Z,V ASCII file

if exist('filevel','var')==0 || isempty(filevel)==1
    [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select velocity model');
else
    pathvel=[];
end
Vfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vp)

aa=load(Vfile);
if size(aa,2)>2
    aa(isnan(aa(:,3)),:)=[];
end

[xnew,ynew]=local2utm(aa(:,1),utm_start,utm_end);

[~,newfile] = fileparts(Vfile);
Vfile_new = [pathvel,newfile,'_UTM','.dat'];
if size(aa,2)>2
    dlmwrite(Vfile_new,[xnew,ynew,aa(:,2),log10(aa(:,3))],'delimiter','\t','precision','%6.2f');
else
    dlmwrite(Vfile_new,[xnew,ynew,log10(aa(:,2))],'delimiter','\t','precision','%6.2f');
end
