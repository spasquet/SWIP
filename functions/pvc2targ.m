function lmaxpick=pvc2targ(pvcstruct,dir_pick,nametarg,wave,sampling,resampvec,flim,maxerr)

%%% S. Pasquet - V22.05.04
% Convert .pvc ASCII file in .target dinver file
% lmaxpick=pvc2targ(pvcstruct,dir_pick,nametarg,wave,sampling,resampvec,flim,maxerr)
% Handles dinver > 2

version = geopsy_version();

if nargin<7
    flim=0;
end
if nargin<8
    maxerr=0.5;
end

if strcmp(wave,'R')==1
    wave='Rayleigh';
elseif strcmp(wave,'L')==1
    wave='Love';
end

% if lininv==1
%     named=[num2str(XmidT(ix),['%10.',num2str(bcscale),'f']),...
%         '.contents.d'];
%     fid2=fopen(named,'w');
% end

nmode=length(pvcstruct);
lmaxpick=zeros(1,nmode);
NSP = [];

fid0=fopen('contents.xml','w');
fprintf(fid0,'%s \n','<Dinver>');
fprintf(fid0,'%s \n','<pluginTag>DispersionCurve</pluginTag>');
fprintf(fid0,'%s \n','<pluginTitle>Surface Wave Inversion</pluginTitle>');
fprintf(fid0,'%s \n','<TargetList>');
if version > 2
    fprintf(fid0,'%s \n','<position>0 0 0</position>');
end
if version == 2
    fprintf(fid0,'%s \n',' <ModalCurveTarget type="dispersion">');
else
    fprintf(fid0,'%s \n',' <DispersionTarget type="dispersion">');
end
fprintf(fid0,'%s \n','<selected>true</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
if version == 2
    fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
else
    fprintf(fid0,'%s \n','<misfitType>L2_LogNormalized</misfitType>');
end
for ip=1:nmode
    pvcfile=pvcstruct(ip).name;
    m=str2double(pvcfile(end-4)); % Mode number
    Vprev=load(fullfile(dir_pick,pvcfile));
    % Resample in lambda or frequency
    Vprev(Vprev(:,1)==0,:)=[];
    [freqresamp,vresamp,deltaresamp]=resampvel(Vprev(:,1),Vprev(:,2),...
        Vprev(:,3),resampvec,sampling);
    vresamp=vresamp(freqresamp>flim);
    deltaresamp=deltaresamp(freqresamp>flim);
    freqresamp=freqresamp(freqresamp>flim);
    if isempty(freqresamp)==0 && length(freqresamp)>1
        lmaxpick(ip)=max(vresamp./freqresamp);
    else
        lmaxpick(ip)=NaN;
    end
    nsp=length(vresamp);
    if nsp <= 1
        continue
    end
    
    NSP=NSP+nsp;
    fprintf(fid0,'%s \n','<ModalCurve>');
    fprintf(fid0,'%s \n',['<name>',pvcfile,'</name>']);
    fprintf(fid0,'%s \n',['<log>',num2str(nsp),...
        ' samples loaded from file ',pvcfile,'</log>']);
    if version > 2
        fprintf(fid0,'%s \n','<enabled>true</enabled>');
    end
    fprintf(fid0,'%s \n','<Mode>');
    fprintf(fid0,'%s \n','<slowness>Phase</slowness>');
    if version == 2
        fprintf(fid0,'%s \n',['<polarisation>',...
            wave,'</polarisation>']);
    else
        fprintf(fid0,'%s \n',['<polarization>',...
            wave,'</polarization>']);
    end
    fprintf(fid0,'%s \n','<ringIndex>0</ringIndex>');
    fprintf(fid0,'%s %i %s \n','<index>',m,'</index>');
    fprintf(fid0,'%s \n','</Mode>');
    for ic=1:nsp
        %         stdev=1/(vresamp(ic)-deltaresamp(ic))-...
        %             1/(vresamp(ic)+deltaresamp(ic));
        
        stdev=deltaresamp(ic)/((vresamp(ic))^2-deltaresamp(ic)^2);
        if stdev>maxerr*(1/vresamp(ic))
            stdev=maxerr*(1/vresamp(ic));
        end
        
        if version == 2
            fprintf(fid0,'%s \n','<StatPoint>');
        else
            fprintf(fid0,'%s \n','<RealStatisticalPoint>');
        end
        fprintf(fid0,'%s %4.18f %s\n','<x>',...
            freqresamp(ic),'</x>');
        fprintf(fid0,'%s %4.18f %s\n',' <mean>',...
            1/vresamp(ic),'</mean>');
        fprintf(fid0,'%s %4.18f %s\n','<stddev>',...
            stdev,'</stddev>');
        fprintf(fid0,'%s \n','<weight>1</weight>');
        fprintf(fid0,'%s \n','<valid>true</valid>');
        if version == 2
            fprintf(fid0,'%s \n',' </StatPoint>');
        else
            fprintf(fid0,'%s \n','</RealStatisticalPoint>');
        end        
        %         if lininv==1
        %             fprintf(fid2,'%s %d %4.6f %4.6f %4.6f\n',...
        %                 'SURF96 R C X',m,1/freqresamp(ic),...
        %                 vresamp(ic)/1000,deltaresamp(ic)/1000);
        %         end
    end
    fprintf(fid0,'%s \n','</ModalCurve>');
end

if all(isnan(lmaxpick))
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   No sample in the range of resampvec - Target file empty');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    fclose(fid0);
    unix_cmd(['rm -rf ',nametarg]);
    return
end

if NSP>100
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n              Warning: Number of samples > 100');
    fprintf('\n        Dinver will ask confirmation for each inversion');
    fprintf('\n   Run module C with verbose = 1 or decrease number of sample');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
end

if version == 2
    fprintf(fid0,'%s \n','</ModalCurveTarget>');
else
    fprintf(fid0,'%s \n','</DispersionTarget>');
end

fprintf(fid0,'%s \n','<AutocorrTarget>');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_NormalizedBySigmaOnly</misfitType>');
fprintf(fid0,'%s \n','<AutocorrCurves>');
fprintf(fid0,'%s \n','</AutocorrCurves>');
fprintf(fid0,'%s \n','</AutocorrTarget>');
fprintf(fid0,'%s \n','<ModalCurveTarget type="ellipticity">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_LogNormalized</misfitType>');
fprintf(fid0,'%s \n','</ModalCurveTarget>');
fprintf(fid0,'%s \n','<ValueTarget type="ellipticity peak">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
fprintf(fid0,'%s \n','<StatValue>');
fprintf(fid0,'%s \n','<mean>0</mean>');
fprintf(fid0,'%s \n','<stddev>0</stddev>');
fprintf(fid0,'%s \n','<weight>1</weight>');
fprintf(fid0,'%s \n','<valid>false</valid>');
fprintf(fid0,'%s \n','</StatValue>');
fprintf(fid0,'%s \n','</ValueTarget>');
fprintf(fid0,'%s \n','<RefractionTarget type="Vp">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
fprintf(fid0,'%s \n','</RefractionTarget>');
fprintf(fid0,'%s \n','<RefractionTarget type="Vs">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
fprintf(fid0,'%s \n','</RefractionTarget>');
if version > 2
    fprintf(fid0,'%s \n','<MagnetoTelluricTarget>');
    fprintf(fid0,'%s \n','<selected>false</selected>');
    fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
    fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
    fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
    fprintf(fid0,'%s \n','</MagnetoTelluricTarget>');
end
fprintf(fid0,'%s \n','</TargetList>');
fprintf(fid0,'%s \n','</Dinver>');

fclose(fid0);
% if lininv==1
%     fclose(fid2);
% end
% Conversion in target

% Conversion in param
if ismac == 1
    unix_cmd(['gtar czf ',nametarg,' contents.xml']);
elseif ispc == 1
    if ispc_wsl == 0
        unix_cmd(['tar czf ',nametarg,' contents.xml --force-local']);
    else
        unix_cmd(['tar -czf ',nametarg,' "contents.xml"']);
    end
else
    unix_cmd(['tar czf ',nametarg,' contents.xml']);
end
delete('contents.xml');

% if ismac == 1
%     unix(['gtar czf ',nametarg,' contents.xml']);
% elseif ispc == 1 && ispc_wsl == 0
%     unix(['tar czf ',nametarg,' contents.xml --force-local']);
% else
%     unix(['tar czf ',nametarg,' contents.xml']);
% end
% delete('contents.xml');

nametarg=strrep(nametarg,'\','\\');
fprintf(['\n  Target file saved as ',nametarg,'\n']);

end
