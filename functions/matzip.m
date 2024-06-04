function matzip(compress,files,method,del,output,verbose)

%%% S. Pasquet - V22.05.04
% GZIP in matlab
% matzip(compress,files,method,del,output,verbose)
wsl = ispc_wsl;

if nargin<5
    output=[];
end

if nargin<6
    verbose=0;
end

if verbose==1
    tic
end
if compress==1
    if method==1
        gzip(files);
        if del==1
            delete(files);
        end
    elseif method==2
        tar(output,files);
        if del==1
            delete(files);
        end
    elseif method==3
        zip(output,files);
        if del==1
            delete(files);
        end
    elseif method==4
        [~,~]=unix_cmd(['7z a -aoa ',output,' ',files],wsl);
        if del==1
            [~,~]=unix_cmd(['rm -rf ',files],wsl);
        end
    elseif method==5
        if del==1
            [~,~]=unix_cmd(['gzip -f ',files],wsl);
        else
            [~,~]=unix_cmd(['gzip -k -f ',files],wsl);
        end
    elseif method==6
        if del==1
            [~,~]=unix_cmd(['bzip2 -f ',files],wsl);
        else
            [~,~]=unix_cmd(['bzip2 -k -f ',files],wsl);
        end
    end
elseif compress==0
    if method==1
        gunzip(files);
        if del==1
            delete(files);
        end
    elseif method==2
        untar(files);
        if del==1
            delete(files);
        end
    elseif method==3
        unzip(files,output);
        if del==1
            delete(files);
        end
    elseif method==4
        [~,~]=unix_cmd(['7z e ',files],wsl);
        if del==1
            [~,~]=unix_cmd(['rm -rf ',files],wsl);
        end
    elseif method==5
        if del==1
            [~,~]=unix_cmd(['gunzip -f ',files],wsl);
        else
            unix_cmd(['gunzip -k -f ',files],wsl);
        end
    elseif method==6
        if del==1
            [~,~]=unix_cmd(['bunzip2 -f ',files],wsl);
        else
            [~,~]=unix_cmd(['bunzip2 -k -f ',files],wsl);
        end
    end
else
    fprintf('\n Wrong command\n');
end
if verbose==1
    toc
end