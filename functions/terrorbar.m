function herrorbars=terrorbar(varargin)
%function herrorbars=terrorbar(x,val,lowererror,uppererror,errorbarwidth,errorbarunits)
%
%=========================
% terrorbar.m 
% Draws error bars (just the error bars, not the lines) whose size can be
% controlled (which was otherwise a challenge in versions 2014b onwards)
%
% Trevor Agus, 31st July 2015
%=========================
%
%
%function herrorbars=terrorbar(x,val,symmetricerrorsize,errorbarwidth,errorbarunits)
%function herrorbars=terrorbar(x,val,symmetricerrorsize,errorbarwidth) %assumes width is in x-axis units
%function herrorbars=terrorbar(x,val,lowererror,uppererror,errorbarwidth)
% 
%%after error lines have been created, these also work to update them:
%function terrorbar(herrorbars,errorbarwidth,errorbarunits)
%function terrorbar(errorbarwidth,errorbarunits)
%function terrorbar(errorbarwidth)
%function terrorbar %e.g. to redraw the bars after resizing the axis
%
% Examples:
% %Basic function: draws error bars
% >> herrorbars=terrorbar([1 2 3],[2 4 6],[.3 .5 .7],[.5 .7 .9],.4,'centi');
% >> herrorbars=terrorbar([1 2 3],[2 4 6],[.5 .7 .9],.4,'centi'); %symmetric
% >> herrorbars=terrorbar([1 2 3],[2 4 6],[.3 .5 .7],[.5 .7 .9],.25,'units'); %proportional to x-axis units
% >> herrorbars=terrorbar([1 2 3],[2 4 6],[.3 .5 .7],[.5 .7 .9],.25); %...which is the default
% >> herrorbars=terrorbar([1 2 3],[2 4 6],[.3 .5 .7],.25); %symmetric
%
% %The error bar handles are simple handles to Lines, so you can easily tweak the properties after:
% >> set(herrorbars,'LineWidth',2,'Color','r')
%
% %You can update the width of the error bars after you've drawn them:
% >> terrorbar(herrorbars,.2,'centi');
% >> terrorbar(herrorbars,.4); %doesn't change units
%
% %This will affect the size of all error bars in the current figure:
% 
% >> terrorbar(.2,'centi')
%
% %Should update as figure size changes, but will only automatically update
% %with axis-size changes if 'units' or 'norm' are used (as opposed to 'centimetres'). 
% %Changing the xlimits can also pose problems for absolute measurements.
% %However, the error bars' width can always be tweaked to the correct widths by simply calling the same command without arguments:
%
% >> terrorbar

global TERRORLOGMEMORY
TERRORLOGMEMORY.axissizechanged=false;

switch(nargin)
    case {4 5 6} %normal initial plotting of errorbars
        herrorbars=internal_ploterrorbars(varargin{:});
    case {0} %this will be called if the figure is resized, or can be called manually if the axis is resized
        internal_tweakallerrorbars
    case {1 2 3} %to update the size later
        if ~isnumeric(varargin{1})
            internal_tweakspecificerrorbars(varargin{:}) %to update specified error bars that have already been created
            herrorbars=varargin{1};
        else
            internal_tweakallerrorbars(varargin{:}) %will update all error bars in a figure; two inputs max
        end;
    otherwise
        error('Unexpected switch')
end;
if TERRORLOGMEMORY.axissizechanged %Oops: XLim must have changed mid-way through tweaking the errorbar widths, potentially leading to errorbars of different widths
    if exist('herrorbars','var')
        terrorbar(herrorbars);
    else
        terrorbar; %running this program (recursively) tweaks the error bar widths again, taking the latest XLim range into account
    end;
end;

function h=internal_ploterrorbars(varargin)

%hard-wired variables
initbarwidth=0.00001; %will tweak this later

%sort inputs
switch(nargin)
    case {6} %all arguments
        x=varargin{1};
        val=varargin{2};
        lowererror=varargin{3};
        uppererror=varargin{4};
        errorbarwidth=varargin{5};
        errorbarunits=varargin{6};
    case {4} %
        x=varargin{1};
        val=varargin{2};
        lowererror=varargin{3};
        uppererror=varargin{3}; %assumes symmetric bars
        errorbarwidth=varargin{4};
        errorbarunits='units'; %defaults to units of the x axis
    case {5}
        if ischar(varargin{5})
            x=varargin{1};
            val=varargin{2};
            lowererror=varargin{3};
            uppererror=varargin{3}; %assumes symmetric bars
            errorbarwidth=varargin{4};
            errorbarunits=varargin{5};
        else
            x=varargin{1};
            val=varargin{2};
            lowererror=varargin{3};
            uppererror=varargin{4};
            errorbarwidth=varargin{5};
            errorbarunits='units'; %defaults to units of the x axis            
        end;
    otherwise
        error('Expected 4, 5, or 6 arguments')
end;

%We need 'hold on' here, but don't want to mess up other settings
oldhold=ishold; %will reset this later

%plot horizontal bars. We'll worry about their correct width in a separate
%subfunction, just to avoid duplication of code
for ii=length(x):-1:1
    h(ii+2*length(x))=plot([x(ii)-initbarwidth x(ii)+initbarwidth],[1 1]*(val(ii)+uppererror(ii)),'k-');
    if ii==length(x); hold on; end;
    h(ii+length(x))=plot([x(ii)-initbarwidth x(ii)+initbarwidth],[1 1]*(val(ii)-lowererror(ii)),'k-');
end;
set(h(length(x)+1:3*length(x)),'UserData',{828399489 errorbarwidth errorbarunits}) %useful for resizing later; 828399489 is just an arbitrary code that should be easy to distinguish from any other UserData in an image

%plot the vertical lines of the error bars
for ii=length(x):-1:1
    h(ii)=plot([1 1]*x(ii),[val(ii)-lowererror(ii) val(ii)+uppererror(ii)],'k-');
end;

%reset whatever the 'hold' was
if oldhold; 
    hold on; 
else
    hold off
end;

%now set the width of the error bars to what they should have been
internal_tweakspecificerrorbars(h,errorbarwidth,errorbarunits)

%the following line means that changing the Figure size will redraw the error bars according to the current rules
set(gcf,'SizeChangedFcn','terrorbar') %calls back main function from this file
%...however if the axis is resized, this program won't know about it unless you run "terrorbars" [no arguments necessary]

function internal_tweakallerrorbars(varargin)
handlestosearch=gcf;
currenthandle=0;
while currenthandle<length(handlestosearch)
    currenthandle=currenthandle+1;
    handlestosearch=[handlestosearch; get(handlestosearch(currenthandle),'Children')]; %#ok<AGROW>
    userdata=get(handlestosearch(currenthandle),'UserData');
    if ~isempty(userdata)&&iscell(userdata)&&userdata{1}==828399489
        if nargin>0
            internal_tweakanerrorbar(handlestosearch(currenthandle),varargin{:});
        else
            internal_tweakanerrorbar(handlestosearch(currenthandle));
        end;
    end;
end;

function internal_tweakspecificerrorbars(h,varargin)
for ii=1:length(h)
    internal_tweakanerrorbar(h(ii),varargin{:})
end;

function internal_tweakanerrorbar(h,errorbarwidth,errorbarunits)
global TERRORLOGMEMORY

xdata=get(h,'XData');
userdata=get(h,'UserData');
if isempty(userdata)||~iscell(userdata)||userdata{1}~=828399489
    %not a horizontal errorbar. Do nothing
    return
end;

%update the line's information, if appropriate; otherwise update 'userdata'
if exist('errorbarwidth','var')
    userdata{2}=errorbarwidth;
else
    errorbarwidth=userdata{2};
end;
if exist('errorbarunits','var')
    userdata{3}=errorbarunits;
else
    errorbarunits=userdata{3};
end;
set(h,'UserData',userdata)

if strcmp(get(gca,'XScale'),'log')
    error('terrorlog.m will behave strangely with a logarithmic x-axis: plot log(x) linearly, tweak the code, or contact trevor@earcatching.com to get him to finish writing terrorlog.m')
end;

centrepoint=mean(xdata);
halfwidth=errorbarwidth/2;
switch(errorbarunits) %errorbarunits
    case {'unit' 'units'}
        xdata=[centrepoint-halfwidth centrepoint+halfwidth];
    case {'pixel' 'pixe' 'pixels'}
        %get current width of current axis in pixels
        oldunits=get(gca,'Units');
        set(gca,'Units','pixels')
        pos=get(gca,'Position');
        set(gca,'Units',oldunits)
        width=pos(3);
        
        %convert specifications to plottable values on the graph
        xlimits=xlim;
        xrange=xlimits(2)-xlimits(1);
        conversationrate=xrange/width;
        xdata=[centrepoint-halfwidth*conversationrate centrepoint+halfwidth*conversationrate];
    case {'centi' 'cent' 'centimetres'}
        %get current width of current axis in centimetres
        oldunits=get(gca,'Units');
        set(gca,'Units','centi')
        pos=get(gca,'Position');
        set(gca,'Units',oldunits)
        width=pos(3);
        
        %convert specifications to plottable values on the graph
        xlimits=xlim;
        xrange=xlimits(2)-xlimits(1);
        conversationrate=xrange/width;
        xdata=[centrepoint-halfwidth*conversationrate centrepoint+halfwidth*conversationrate];
    case {'inches' 'inch'}
        %get current width of current axis in inches
        oldunits=get(gca,'Units');
        set(gca,'Units','inches')
        pos=get(gca,'Position');
        set(gca,'Units',oldunits)
        width=pos(3);
        
        %convert specifications to plottable values on the graph
        xlimits=xlim;
        xrange=xlimits(2)-xlimits(1);
        conversationrate=xrange/width;
        xdata=[centrepoint-halfwidth*conversationrate centrepoint+halfwidth*conversationrate];
    case {'points' 'point' 'poin'}
        %get current width of current axis in points
        oldunits=get(gca,'Units');
        set(gca,'Units','points')
        pos=get(gca,'Position');
        set(gca,'Units',oldunits)
        width=pos(3);
        
        %convert specifications to plottable values on the graph
        xlimits=xlim;
        xrange=xlimits(2)-xlimits(1);
        conversationrate=xrange/width;
        xdata=[centrepoint-halfwidth*conversationrate centrepoint+halfwidth*conversationrate];
    case {'norm' 'normalised' 'normalized'} %with respect to current axis
        %convert specifications to plottable values on the graph
        xlimits=xlim;
        xrange=xlimits(2)-xlimits(1);
        conversationrate=xrange; %/1 where "1" is the width of the axis
        xdata=[centrepoint-halfwidth*conversationrate centrepoint+halfwidth*conversationrate];
    otherwise
        userdata{3} %#ok<NOPRT>
        error('Unfamiliar units')
end;

%Tweak width
previousrange=get(gca,'XLim');
set(h,'XData',xdata)
if any(previousrange~=get(gca,'XLim')) %XLim can be automatically changed when width of errorbars are changed
    TERRORLOGMEMORY.axissizechanged=true;
end;