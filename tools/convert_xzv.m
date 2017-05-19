function convert_xzv(filevel,xmin,xmax,ymin,ymax)

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

% keyboard

x=aa(:,1);
y=zeros(size(x));

% figure
% plot(x,y,'k+');
% axis equal

teta=atan((ymax-ymin)/(xmax-xmin));

xnew=x*cos(teta)+y*sin(teta)+xmin;
ynew=y*cos(teta)+x*sin(teta)+ymin;

% figure
% plot(xnew,ynew,'kx')
% axis equal

[~,newfile] = fileparts(Vfile);
Vfile_new = [pathvel,newfile,'_UTM','.dat'];
if size(aa,2)>2
    dlmwrite(Vfile_new,[xnew,ynew,aa(:,2),aa(:,3)],'delimiter','\t','precision','%6.2f');
else
    dlmwrite(Vfile_new,[xnew,ynew,aa(:,2)],'delimiter','\t','precision','%6.2f');
end
    
% 
% 	xA = 544446.587;
% 	yA = 4939636.387;
%         
% 	xB = 544408.884;
% 	yB = 4939644.172;
%     
%     xC = 544403.157;
%     yC = 4939624.181;
%     
%     xD = 544440.798;
%     yD = 4939616.171;
%     
%     plot([xA,xC,xB,xD],[yA,yC,yB,yD],'kx-')
%     hold on
%     plot((xC+xA)/2,(yC+yA)/2,'b+')
%     
%     xtrans = (xC+xA)/2;
%     ytrans = (yC+yA)/2;
%     
%     teta=180+180*atan((yB-yA)/(xB-xA))/pi;
