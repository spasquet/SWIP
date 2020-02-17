function plot_seismo(seismofile,tMIN,tMAX,imgform,imgres,fs,fldr,scal,clip,perc,xticks,tticks,xMIN,xMAX,imgsize)

%%% S. Pasquet - V16.11.22
% Quick plot of seismogram
% plot_seismo(seismofile,tMIN,tMAX,imgform,imgres,fs,fldr,scal,clip,perc,xticks,tticks,xMIN,xMAX,imgsize)

if exist('seismofile','var')==0 || isempty(seismofile)==1
    [seismofile,seismopath]=uigetfile('*.su','Select seismogram file');
    if seismofile==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a seismogram file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
else
    seismopath=[];
end

run('SWIP_defaultsettings')
% keyboard
if exist('fldr','var')==0 || isempty(fldr)==1
    % Get fldr in case two shots at the same position
    com1=sprintf('sugethw < %s key=fldr output=geom | uniq',fullfile(seismopath,seismofile));
    [~,fldr]=unix(com1);
    fldr=str2num(fldr);
end
if exist('scal','var')==0 || isempty(scal)==1
    scal=1;
end
if exist('clip','var')==0 || isempty(clip)==1
    clip=1;
end
if exist('perc','var')==0 || isempty(perc)==1
    perc=97;
end

if exist('imgsize','var')==0 || isempty(imgsize)==1
    imgsize=[0 0 25 15];
end

[~,xsca]=unix(['sugethw < ',fullfile(seismopath,seismofile),' key=scalco output=geom | uniq']);
xsca=abs(str2double(xsca));
if xsca == 0
   xsca = 1; 
end

tMIN = -20;
tMAX = 200;

for i=1:length(fldr)
    com1=sprintf('suwind < %s key=fldr min=%d max=%d > tmp.su',fullfile(seismopath,seismofile),fldr(i),fldr(i));
    unix(com1);
    [seismomat,tseis,xseis]=seismo2dat(fullfile(seismopath,'tmp.su'),0);
    seismomat_norm = bsxfun(@rdivide,seismomat,max(abs(seismomat),[],2));
    delete(fullfile(seismopath,'tmp.su'));
    fig=plot_wiggle(showplot,-seismomat_norm',xseis/xsca,tseis*1000,scal,clip,perc,...
        fs,'Gx (m)','Time (ms)',[xMIN xMAX],[tMIN tMAX],xticks,tticks,imgsize,[]);
    file1=[num2str(fldr(i)),'.seismo.',imgform];
    save_fig(fig,file1,imgform,imgres,1);
end

end