function paramname=mod2param(nlay,nsublay,thmin,thmax,shape,lvz,...
    Vpmin,Vpmax,Numin,Numax,Vsmin,Vsmax,Rhomin,Rhomax,...
    Vplink,Nulink,Vslink,Rholink,paramname)
%
% %% S. Pasquet - V16.9.13
% Create parameterization for dinver
% nlay = Nb of layers (including half-space)
%
% nsublay = Nb of sublayers per layer
% shape = Shape of the velocity variation with depth
% (1='Uniform', 2='Linear', 3='LinearIncrease', 4='LinearDecrease', 
% 5='PowerLaw')
%
% thmin = Minimum thickness (m)
% thmax = Maximum thickness (m)
% Vpmin = Minimum Vp (m/s)
% Vpmax = Maximum Vp (m/s)
% Numin = Minimum Poisson's ratio
% Numax = Maximum Poisson's ratio
% Vsmin = Minimum Vs (m/s)
% Vsmax = Maximum Vs (m/s)
% Rhomin = Minimum Density (kg/m3)
% Rhomax = Maximum Density (kg/m3)
%
% lvz = authorize low velocity zone (=1) or not (=0)
%
% Vplink = Vp linked to Vs (0 = not linked, 1 = linked)
% Nulink = Nu linked to Vs (0 = not linked, 1 = linked)
% Vslink = Vs linked to Vs (0 = not linked, 1 = linked)
% Rholink = Rho linked to Vs (0 = not linked, 1 = linked)
%
% Input other than nlay and paramname can be either scalar or vector with 
% the length of nlay
%
% paramname = Specific filename (if none entered, auto filename)

% Repeat values for nlay if only one value entered
if length(nsublay)==1
    nsublay=repmat(nsublay,nlay,1);
elseif length(nsublay)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(nsublay)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(shape)==1
    shape=repmat(shape,nlay,1);
elseif length(shape)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(shape)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(thmin)==1
    thmin=repmat(thmin,nlay,1);
elseif length(thmin)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(thmin)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(thmax)==1
    thmax=repmat(thmax,nlay,1);
elseif length(thmax)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(thmax)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Vpmin)==1
    Vpmin=repmat(Vpmin,nlay,1);
elseif length(Vpmin)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Vpmin)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Vpmax)==1
    Vpmax=repmat(Vpmax,nlay,1);
elseif length(Vpmax)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Vpmax)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Numin)==1
    Numin=repmat(Numin,nlay,1);
elseif length(Numin)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Numin)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Numax)==1
    Numax=repmat(Numax,nlay,1);
elseif length(Numax)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Numax)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Vsmin)==1
    Vsmin=repmat(Vsmin,nlay,1);
elseif length(Vsmin)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Vsmin)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Vsmax)==1
    Vsmax=repmat(Vsmax,nlay,1);
elseif length(Vsmax)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Vsmax)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Rhomin)==1
    Rhomin=repmat(Rhomin,nlay,1);
elseif length(Rhomin)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Rhomin)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Rhomax)==1
    Rhomax=repmat(Rhomax,nlay,1);
elseif length(Rhomax)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Rhomax)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(lvz)==1
    lvz=repmat(lvz,nlay,1);
elseif length(lvz)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(lvz)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Vplink)==1
    Vplink=repmat(Vplink,nlay,1);
elseif length(Vplink)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Vplink)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Nulink)==1
    Nulink=repmat(Nulink,nlay,1);
elseif length(Nulink)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Nulink)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Vslink)==1
    Vslink=repmat(Vslink,nlay,1);
elseif length(Vslink)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Vslink)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end
if length(Rholink)==1
    Rholink=repmat(Rholink,nlay,1);
elseif length(Rholink)~=nlay
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   Error: length(Rholink)~=nlay');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    paramname=[];
    return
end

% Parameter structure
param(1).min=Vpmin;
param(2).min=Numin;
param(3).min=Vsmin;
param(4).min=Rhomin;
param(1).max=Vpmax;
param(2).max=Numax;
param(3).max=Vsmax;
param(4).max=Rhomax;
param(2).link=Nulink;
param(3).link=Vslink;
param(4).link=Rholink;
param(1).link=Vplink;
param(2).link=Nulink;
param(3).link=Vslink;
param(4).link=Rholink;
param(1).long='Compression-wave velocity';
param(2).long='Poisson''s Ratio';
param(3).long='Shear-wave velocity';
param(4).long='Density';
param(1).short='Vp';
param(2).short='Nu';
param(3).short='Vs';
param(4).short='Rho';

for i=1:nlay
    if shape(i)==1
        layer(i).shape='Uniform';
    elseif shape(i)==2
        layer(i).shape='Linear';
    elseif shape(i)==3
        layer(i).shape='LinearIncrease';
    elseif shape(i)==4
        layer(i).shape='LinearDecrease';
    elseif shape(i)==5
        layer(i).shape='PowerLaw';
    end
    if lvz(i)==1
        layer(i).lvz='false';
    else
        layer(i).lvz='true';
    end
    layer(i).thmin=thmin(i);
    layer(i).thmax=thmax(i);
    layer(i).nsublay=nsublay(i);
end

fid0=fopen('contents.xml','w');
fprintf(fid0,'%s \n','<Dinver>');
fprintf(fid0,'%s \n','<pluginTag>DispersionCurve</pluginTag>');
fprintf(fid0,'%s \n','<pluginTitle>Surface Wave Inversion</pluginTitle>');
fprintf(fid0,'%s \n','<ParamGroundModel>');

for j=1:4
    fprintf(fid0,'%s \n','<ParamProfile>');
    fprintf(fid0,'%s \n','<type>Param</type>');
    fprintf(fid0,'%s \n',['<longName>',param(j).long,'</longName>']);
    fprintf(fid0,'%s \n',['<shortName>',param(j).short,'</shortName>']);
    if j==1
        fprintf(fid0,'%s \n',['<unit>','m/s','</unit>']);
        fprintf(fid0,'%s \n','<defaultMinimum>200</defaultMinimum>');
        fprintf(fid0,'%s \n','<defaultMaximum>5000</defaultMaximum>');
        fprintf(fid0,'%s \n','<defaultCondition>LessThan</defaultCondition>');
    elseif j==2
        fprintf(fid0,'%s \n','<unit/>');
        fprintf(fid0,'%s \n','<defaultMinimum>0.1</defaultMinimum>');
        fprintf(fid0,'%s \n','<defaultMaximum>0.5</defaultMaximum>');
        fprintf(fid0,'%s \n','<defaultCondition>LessThan</defaultCondition>');
    elseif j==3
        fprintf(fid0,'%s \n',['<unit>','m/s','</unit>']);
        fprintf(fid0,'%s \n','<defaultMinimum>50</defaultMinimum>');
        fprintf(fid0,'%s \n','<defaultMaximum>2500</defaultMaximum>');
        fprintf(fid0,'%s \n','<defaultCondition>LessThan</defaultCondition>');
    elseif j==4
        fprintf(fid0,'%s \n',['<unit>','kg/m3','</unit>']);
        fprintf(fid0,'%s \n','<defaultMinimum>1800</defaultMinimum>');
        fprintf(fid0,'%s \n','<defaultMaximum>1800</defaultMaximum>');
        fprintf(fid0,'%s \n','<defaultCondition>LessThan</defaultCondition>');
    end
    for i=1:nlay
        fprintf(fid0,'%s \n',['<ParamLayer name="',param(j).short,num2str(i-1),'">']);
        if i<nlay
            fprintf(fid0,'%s \n',['<shape>',layer(i).shape,'</shape>']);
        elseif i==nlay
            fprintf(fid0,'%s \n',['<shape>','Uniform','</shape>']);
        end
        if i==1
            fprintf(fid0,'%s \n',...
                ['<lastParamCondition>','true','</lastParamCondition>']);
        else
            if j==2
                fprintf(fid0,'%s \n',...
                    ['<lastParamCondition>','false','</lastParamCondition>']);
            else
                fprintf(fid0,'%s \n',...
                    ['<lastParamCondition>',layer(i).lvz,'</lastParamCondition>']);
            end
        end
        fprintf(fid0,'%s \n',['<nSubayers>',num2str(layer(i).nsublay),'</nSubayers>']);
        fprintf(fid0,'%s \n',['<topMin>',num2str(param(j).min(i)),'</topMin>']);
        fprintf(fid0,'%s \n',['<topMax>',num2str(param(j).max(i)),'</topMax>']);
        if strcmp(layer(i).shape,'Linear')==1 || strcmp(layer(i).shape,'LinearIncrease')==1 ...
                || strcmp(layer(i).shape,'LinearDecrease')==1 || strcmp(layer(i).shape,'PowerLaw')==1
            fprintf(fid0,'%s \n',['<bottomMin>',num2str(param(j).min(i)),'</bottomMin>']);
            fprintf(fid0,'%s \n',['<bottomMax>',num2str(param(j).max(i)),'</bottomMax>']);
        end
        if param(j).link(i)==0
            fprintf(fid0,'%s \n',['<linkedTo>','Not linked','</linkedTo>']);
        else
            fprintf(fid0,'%s \n',['<linkedTo>Vs',num2str(i-1),'</linkedTo>']);
        end
        
        if i==1 || i==nlay
            fprintf(fid0,'%s \n',['<isDepth>','true','</isDepth>']);
        else
            fprintf(fid0,'%s \n',['<isDepth>','false','</isDepth>']);
        end
        fprintf(fid0,'%s \n',['<dhMin>',num2str(layer(i).thmin),'</dhMin>']);
        fprintf(fid0,'%s \n',['<dhMax>',num2str(layer(i).thmax),'</dhMax>']);
        fprintf(fid0,'%s \n','</ParamLayer>');
    end
    fprintf(fid0,'%s \n','</ParamProfile>');
end
fprintf(fid0,'%s \n','<ParamSpaceScript>');
fprintf(fid0,'%s \n','<text>');
fprintf(fid0,'%s \n','// model description');
fprintf(fid0,'%s \n','</text>');
fprintf(fid0,'%s \n','</ParamSpaceScript>');
fprintf(fid0,'%s \n','</ParamGroundModel>');
fprintf(fid0,'%s \n','</Dinver>');
fclose(fid0);

if exist('paramname','var')~=1 || isempty(paramname)==1
    paramname=['C',num2str(nlay)];
    for i=1:nlay
        if i==1 || length(unique(layer(i).thmin))~=1 && length(unique(layer(i).thmax))~=1
            if strcmp(layer(i).shape,'Uniform')==1
                paramname=[paramname,'.U.Th',num2str(layer(i).thmin),'_',...
                    num2str(layer(i).thmax)];
            elseif strcmp(layer(i).shape,'Linear')==1
                paramname=[paramname,'.L.Th',num2str(layer(i).thmin),'_',...
                    num2str(layer(i).thmax),'ns',num2str(layer(i).nsublay)];
            elseif strcmp(layer(i).shape,'LinearIncrease')==1
                paramname=[paramname,'.LI.Th',num2str(layer(i).thmin),'_',...
                    num2str(layer(i).thmax),'ns',num2str(layer(i).nsublay)];
            elseif strcmp(layer(i).shape,'LinearDecrease')==1
                paramname=[paramname,'.LD.Th',num2str(layer(i).thmin),'_',...
                    num2str(layer(i).thmax),'ns',num2str(layer(i).nsublay)];
            elseif strcmp(layer(i).shape,'PowerLaw')==1
                paramname=[paramname,'.PL.Th',num2str(layer(i).thmin),'_',...
                    num2str(layer(i).thmax),'ns',num2str(layer(i).nsublay)];
            end
            for j=1:4
                if param(j).link(i)==0
                    if strcmp(layer(i).lvz,'true')==1
                        paramname=[paramname,'.',param(j).short,num2str(param(j).min(i)),...
                            '_',num2str(param(j).max(i))];
                    else
                        paramname=[paramname,'.',param(j).short,num2str(param(j).max(i)),...
                            '_',num2str(param(j).min(i))];
                    end
                else
                    if strcmp(layer(i).lvz,'true')==1
                        paramname=[paramname,'.',param(j).short,'lk',num2str(param(j).min(i)),...
                            '_',num2str(param(j).max(i))];
                    else
                        paramname=[paramname,'.',param(j).short,'lk',num2str(param(j).max(i)),...
                            '_',num2str(param(j).min(i))];
                    end
                end
            end
        end
    end
    paramname=[paramname,'.param'];
else
    [~,~,extparam]=fileparts(paramname);
    if strcmp(extparam,'.param')==0
        paramname=[paramname,'.param'];
    end
end
% Conversion in param
if ismac == 1
    unix(['gtar czf ',paramname,' contents.xml']);
else
    unix(['tar czf ',paramname,' contents.xml']);
end
delete('contents.xml');


end