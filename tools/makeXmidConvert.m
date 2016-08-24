clear all;
dir_pick=uigetdir('./','Select folder containing pvc files to convert');
pvcstruct=dir(fullfile(dir_pick,'*.pvc'));

for i=1:length(pvcstruct) 
pvcfile=pvcstruct(i).name;
XmidTest(i)=str2double(pvcfile(1:end-7));
end
XmidT=unique(XmidTest);
nW=72; % moyenne si plusieurs (type Ezersky)
dx=0.25;
incerr=1; % nWfac
maxerr=0.5; % ratio
minvelerr=10; % m/s
convert_pvc(XmidT,nW,dx,incerr,maxerr,minvelerr)
