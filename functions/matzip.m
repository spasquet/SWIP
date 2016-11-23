function matzip(compress,files,method,del,output,verbose)

%%% S. Pasquet - V16.11.18
% GZIP in matlab
% matzip(compress,files,method,del,output,verbose)
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
        [~,~]=unix(['7z a -aoa ',output,' ',files]);
        if del==1
            [~,~]=unix(['rm -rf ',files]);
        end
    elseif method==5
        if del==1
            [~,~]=unix(['gzip -f ',files]);
        else
            [~,~]=unix(['gzip -k -f ',files]);
        end
    elseif method==6
        if del==1
            [~,~]=unix(['bzip2 -f ',files]);
        else
            [~,~]=unix(['bzip2 -k -f ',files]);
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
        [~,~]=unix(['7z e ',files]);
        if del==1
            [~,~]=unix(['rm -rf ',files]);
        end
    elseif method==5
        if del==1
            [~,~]=unix(['gunzip -f ',files]);
        else
            unix(['gunzip -k -f ',files]);
        end
    elseif method==6
        if del==1
            [~,~]=unix(['bunzip2 -f ',files]);
        else
            [~,~]=unix(['bunzip2 -k -f ',files]);
        end
    end
else
    fprintf('\n Wrong command\n');
end
if verbose==1
    toc
end