function convert_pvc(XmidT,nW,dx,incerr,maxerr,minvelerr)
%%% S. Pasquet - V16.2.17
% Convert old pvc files (2 columns) in new pvc files (3 columns)

dir_pick=uigetdir('./','Select folder containing pvc files to convert');
pvcstruct=dir(fullfile(dir_pick,'*.pvc'));
dir_pick_new=fullfile(dir_pick,'file.pvc');
mkdir(dir_pick_new);

for i=1:length(pvcstruct)
pvcfile=pvcstruct(i).name;
Xmid=str2double(pvcfile(1:end-7));
Xmidok=XmidT(abs(XmidT-Xmid)<dx*0.5);
pvcfilenew=[num2str(Xmidok),pvcfile(end-6:end)];
pvc2pvc(fullfile(dir_pick,pvcfile),nW,dx,incerr,maxerr,...
    minvelerr,fullfile(dir_pick_new,pvcfilenew))
end
