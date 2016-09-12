function [Xi,Zi,Vi]=plot_xzv(filevel,mask,logscale,ztitle,isoline,blocky,vertex,imgform,imgres,vslim)
%%% S. Pasquet - V16.9.12
% Quick plot of 2D sections

if exist('mask','var')==0 || isempty(mask)==1
    mask=0;
end
if exist('logscale','var')==0 || isempty(logscale)==1
    logscale=0;
end
if exist('ztitle','var')==0 || isempty(ztitle)==1
    ztitle='Velocity (m/s)';
end
if exist('isoline','var')==0 || isempty(isoline)==1
    isoline=[];
end
if exist('vslim','var')==1 && isempty(vslim)==0
    vsMIN=vslim(1); vsMAX=vslim(2);
end
if exist('filevel','var')==0 || isempty(filevel)==1
    [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv'},'Select velocity model');
else
    pathvel=[];
end
Vfile=fullfile(pathvel,filevel); % File with velocity (3 columns X,Z,Vp)
if pathvel==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   No file selected - Abort');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    return
else
    try
        [Vi,Xi,Zi]=readtomo(Vfile,mask); % Read Vp tomo file
    catch
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Invalid model file - Abort');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        return
    end
end

run('SWIP_defaultsettings')

if exist('xMIN','var')==0 || isempty(xMIN)==1 || exist('xMAX','var')==0 || isempty(xMAX)==1
    xMIN=min(min(Xi)); xMAX=max(max(Xi));
end

if exist('zMIN','var')==0 || isempty(zMIN)==1 || exist('zMAX','var')==0 || isempty(zMAX)==1
    zMIN=min(min(Zi)); zMAX=max(max(Zi));
end

if logscale==0
    fig=plot_img([],Xi,Zi,Vi,map7,axetop,0,1,fs,'X (m)',...
        'Altitude (m)',ztitle,[xMIN xMAX],[zMIN zMAX],...
        [vsMIN vsMAX],xticks,zticks,vsticks,[],[],isoline,[25 16 24 12],[],vertex,blocky);
else
    fig=plot_img_log([],Xi,Zi,Vi,map7,axetop,0,1,fs,'X (m)',...
        'Altitude (m)',ztitle,[xMIN xMAX],[zMIN zMAX],...
        [vsMIN vsMAX],xticks,zticks,vsticks,[],[],isoline,[25 16 24 12],[],vertex,blocky);
end

file1=[filevel(1:end-3),imgform];
save_fig(fig,file1,imgform,imgres,1);