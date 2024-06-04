function save_fig(h,filename,format,res,crop,verbose,autosize,transparent)

%%% S. Pasquet - V22.05.04
% Save figure in various format
% save_fig(h,filename,format,res,crop,verbose,autosize)

if exist('verbose','var')==0 || isempty(verbose)==1
    verbose=1;
end

if exist('autosize','var')==0 || isempty(autosize)==1
    autosize=1;
end

if exist('transparent','var')==0 || isempty(transparent)==1
    transparent = 0;
end

if transparent == 0
    rend = '-zbuffer';
else
    rend = '-opengl';
end

matrelease=version('-release');
if str2double(matrelease(1:4))>2014
    h=get(h,'Number');
    if autosize==1
        set(h,'paperpositionmode','auto');
    end
    if strcmpi(format,'fig')==1
        savefig(h,filename);
    elseif strcmpi(format,'pdf')==1
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-painters',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
    else
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-opengl',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
    end
else
    if autosize==1
        set(h,'paperpositionmode','auto');
    end
    if strcmpi(format,'fig')==1
        saveas(h,filename);
    elseif strcmpi(format,'pdf')==1
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-painters',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
    else
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            rend,[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
    end
end

filename_unix = unix_wsl_path(filename);

if crop==1 && strcmpi(format,'fig')~=1
    if strcmpi(format,'pdf')==1
        [test,~]=unix_cmd(['pdfcrop -margins 10 ',filename_unix,' ',filename_unix]);
    else
        if isunix == 1 || ispc_wsl == 1
            [test,~]=unix_cmd(['convert ',filename_unix,' -trim ',filename_unix]);
            [~,~]=unix_cmd(['convert ',filename_unix,' -bordercolor White -border 50 ',filename_unix]);
        else
            [test,~]=unix_cmd(['img_convert ',filename_unix,' -trim ',filename_unix]);
            [~,~]=unix_cmd(['img_convert ',filename_unix,' -bordercolor White -border 50 ',filename_unix]);
        end
    end
end
if verbose==1
    filename=strrep(filename,'\','\\');
    fprintf(['\n      Saved as ',filename,'\n']);
end
end
