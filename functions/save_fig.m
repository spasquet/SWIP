function save_fig(h,filename,format,res,crop,verbose,autosize)

%%% S. Pasquet - V16.6.28
% Save figure in various format

if exist('verbose','var')==0 || isempty(verbose)==1
    verbose=1;
end

if exist('autosize','var')==0 || isempty(autosize)==1
    autosize=1;
end

matrelease=version('-release');
if str2double(matrelease(1:4))>2014
    h=get(h,'Number');
    if autosize==1
        set(h,'paperpositionmode','auto');
    end
    if strcmp(format,'fig')==1
        savefig(h,filename);
    elseif strcmp(format,'pdf')==1
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-painters',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
%         unix(['mv ',[filename,'.tmp'],' -- ',filename]);
    else
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-opengl',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
%         unix(['mv ',[filename,'.tmp'],' -- ',filename]);
    end
else
     if autosize==1
        set(h,'paperpositionmode','auto');
    end
    if strcmp(format,'fig')==1
        saveas(h,filename);
    elseif strcmp(format,'pdf')==1
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-painters',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
%         unix(['mv ',[filename,'.tmp'],' -- ',filename]);
    else
        print(['-f' num2str(h)],['-r',num2str(res)],['-d',format],...
            '-zbuffer',[filename,'.tmp']);
        movefile([filename,'.tmp'],filename);
%         unix(['mv ',[filename,'.tmp'],' -- ',filename]);
    end
end

if crop==1 && strcmp(format,'fig')~=1
    if strcmp(format,'pdf')==1
        [test,~]=unix(['pdfcrop -margins 10 ',filename,' ',filename]);
%         if test~=0
%             fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
%             fprintf('\n   Install pdfcrop to crop PDF files');
%             fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
%         end
    else
        if isunix==1
            [test,~]=unix(['convert ',filename,' -trim ',filename]);
            [~,~]=unix(['convert ',filename,' -bordercolor White -border 50 ',filename]);
        else
            [test,~]=unix(['img_convert ',filename,' -trim ',filename]);
            [~,~]=unix(['img_convert ',filename,' -bordercolor White -border 50 ',filename]);
        end
%         if test~=0
%             fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
%             fprintf(['\n   Install Imagemagick to crop ',format,' files']);
%             fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
%         end
    end
end
if verbose==1
    filename=strrep(filename,'\','\\');
    fprintf(['\n      Saved as ',filename,'\n']);
end
end
