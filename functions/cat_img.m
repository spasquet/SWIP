function cat_img(infiles,imgformat,columns,align,outfile,verbose)

%%% S. Pasquet - V16.11.18
% Concatenate images files in various format
%
% cat_img(infiles,imgformat,columns,align,outfile,verbose)

if exist('infiles','var')==0 || isempty(infiles)==1
    infiles=[];
    [file,pathfile]=uigetfile({'*.png;*.jpg;*.jpeg;*.tiff;*.pdf'},'Select image file to concatenate','multiselect','on');
    for i=1:length(file)
        infiles=[infiles,' ',pathfile,file{i}];
    end
else
    pathfile=[];
end

if pathfile==0
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!');
    fprintf('\n   No file selected - Abort');
    fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    return
end

infiles_cell = str_split(infiles,' ');
nfiles = 0;
for i=1:length(infiles_cell)
    test = dir(infiles_cell{i});
    nfiles = nfiles + length(test);
end


if exist('imgformat','var')==0 || isempty(imgformat)==1
    [~,~,ext]=fileparts(infiles);
    imgformat=ext(2:end);
    if strcmp(imgformat,'pdf')==0 && strcmp(imgformat,'png')==0 && strcmp(imgformat,'tiff')==0 &&...
            strcmp(imgformat,'jpg')==0 && strcmp(imgformat,'jpeg')==0
        imgformat='pdf';
    end
end

if exist('columns','var')==0 || isempty(columns)==1 || columns==0
    columns=1;
end

if exist('outfile','var')==0 || isempty(outfile)==1
    outfile=fullfile(pwd,strcat('outfile.',imgformat));
end

if exist('align','var')==0 || isempty(align)==1
    align='center';
end

if exist('verbose','var')==0 || isempty(verbose)==1
    verbose=1;
end

if strcmp(imgformat,'png')==1 || strcmp(imgformat,'jpg')==1 ||...
        strcmp(imgformat,'jpeg')==1 || strcmp(imgformat,'tiff')==1
    x=['-tile ',num2str(columns),'x'];
    com1=sprintf('montage %s -gravity %s -mode concatenate %s %s',infiles,align,x,outfile);
    if isunix==0
        com1=strrep(com1,'\','/');
    end
    unix(com1);
elseif strcmp(imgformat,'pdf')==1
    x=[num2str(columns),'x',num2str(nfiles-columns+1)];
    %     com1=sprintf('pdfjam -q --noautoscale true --papersize ''{%dcm,%dcm}'' --nup  %s %s -o %s',200,200,x,infiles,outfile);
    com1=sprintf('pdfjam -q --noautoscale true --a0paper --fitpaper true --nup  %s %s -o %s',x,infiles,outfile);
    [~,~]=unix(com1);
    [~,~]=unix(['pdfcrop -margins 10 ',outfile,' ',outfile]);
end
if verbose==1
    outfile=strrep(outfile,'\','\\');
    fprintf(['\n      Saved as ',outfile,'\n']);
end
end