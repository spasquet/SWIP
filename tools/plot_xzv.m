function [Xi,Zi,Vi]=plot_xzv(filevel,mask,logscale,xLIM,xticks,zLIM,zticks,vLIM,...
    vsticks,vtitle,isoline,blocky,vertex,imgform,imgres,ISOvel2,filevel2)

%%% S. Pasquet - V16.11.22
% Quick plot of 2D models from X,Z,V ASCII file
% plot_xzv(filevel,mask,logscale,xLIM,xticks,zLIM,zticks,vLIM,...
%     vsticks,vtitle,isoline,blocky,vertex,imgform,imgres,ISOvel2,filevel2)


if exist('mask','var')==0 || isempty(mask)==1
    mask=0;
end
if exist('logscale','var')==0 || isempty(logscale)==1
    logscale=0;
end
if exist('xLIM','var')==1 && isempty(xLIM)==0
    xMIN=xLIM(1); xMAX=xLIM(2);
end
if exist('zLIM','var')==1 && isempty(zLIM)==0
    zMIN=zLIM(1); zMAX=zLIM(2);
end
if exist('vLIM','var')==1 && isempty(vLIM)==0
    vsMIN=vLIM(1); vsMAX=vLIM(2);
end
if exist('vtitle','var')==0 || isempty(vtitle)==1
    vtitle='Velocity (m/s)';
end
if exist('isoline','var')==0 || isempty(isoline)==1
    isoline=[];
end
if exist('filevel','var')==0 || isempty(filevel)==1
    [filevel,pathvel]=uigetfile({'*.model;*.dat;*.xzv;*.txt'},'Select velocity model');
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
% keyboard
if exist('ISOvel2','var')==1 && isempty(ISOvel2)==0
    if exist('filevel2','var')==0 || isempty(filevel2)==1
        [filevel2,pathvel2]=uigetfile({'*.model;*.dat;*.xzv'},'Select velocity model');
    else
        pathvel2=[];
    end
    Vfile2=fullfile(pathvel2,filevel2); % File with velocity (3 columns X,Z,Vp)
    if pathvel2==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No file selected - Abort');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        return
    else
        try
            [Vi2,Xi2,Zi2]=readtomo(Vfile2,mask); % Read Vp tomo file
        catch
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid model file - Abort');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            return
        end
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
        'Altitude (m)',vtitle,[xMIN xMAX],[zMIN zMAX],...
        [vsMIN vsMAX],xticks,zticks,vsticks,[],[],isoline,[25 16 24 12],[],vertex,blocky);
else
    fig=plot_img_log([],Xi,Zi,Vi,map7,axetop,0,1,fs,'X (m)',...
        'Altitude (m)',vtitle,[xMIN xMAX],[zMIN zMAX],...
        [vsMIN vsMAX],xticks,zticks,vsticks,[],[],isoline,[25 16 24 12],[],vertex,blocky);
end

%%%
if exist('ISOvel2','var')==1 && isempty(ISOvel2)~=1
    hold on;
    if length(ISOvel2)==1
        isoline=[ISOvel2 ISOvel2];
    else
        isoline=ISOvel2;
    end
    [cs,hc]=contour(Xi2,Zi2,Vi2,isoline,'color',[0 0 0],'linewidth',1);
    clabel(cs, hc,'manual','Color', 'k', 'Rotation', 0,'fontsize',12);
    hold off;
end
%%%

file1=[filevel(1:end-3),imgform];
save_fig(fig,file1,imgform,imgres,1);