function a2 = magnify2(f1,zoom_lvl,frame_size)

if exist('f1','var')==0 || isempty(f1)==1
    f1 = gcf;
end
if exist('zoom_lvl','var')==0 || isempty(zoom_lvl)==1
    zoom_lvl = 2;
end
if exist('frame_size','var')==0 || isempty(frame_size)==1
    frame_size = 0.2;
end

set(f1,'WindowButtonMotionFcn', @ButtonMotionCallback);

a1 = get(f1,'CurrentAxes');
a2 = copyobj(a1,f1);

set(f1, ...
    'UserData',[f1,a1,a2], ...
    'Pointer','fullcrosshair', ...
    'CurrentAxes',a2);
set(a2, ...
    'UserData',[zoom_lvl,frame_size], ...  %magnification, frame size
    'Color',get(a1,'Color'), ...
    'Box','on');
xlabel(''); ylabel(''); zlabel(''); title('');
set(get(a2,'Children'), ...
    'LineWidth', 2);
set(a1, ...
    'Color',get(a1,'Color')*0.95);
set(f1, ...
    'CurrentAxes',a1);
ButtonMotionCallback(f1);

function ButtonMotionCallback(src,eventdata)
   H = get(src,'UserData');
   if ~isempty(H)
      f1 = H(1); a1 = H(2); a2 = H(3);
      a2_param = get(a2,'UserData');
      f_pos = get(f1,'Position');
      a1_pos = get(a1,'Position');

      [f_cp, a1_cp] = pointer2d(f1,a1);
      set(a2,'Position',[(f_cp./f_pos(3:4)) 0 0]+a2_param(2)*a1_pos(3)*[-0.75 -1.5 1.5 2.25]); %% Figure size
      a2_pos = get(a2,'Position');

   	set(a2,'XLim',a1_cp(1)+(1/a2_param(1))*(a2_pos(3)/a1_pos(3))*diff(get(a1,'XLim'))*[-0.33 0.33]); % Figure axis
   	set(a2,'YLim',a1_cp(2)+(1/a2_param(1))*(a2_pos(4)/a1_pos(4))*diff(get(a1,'YLim'))*[-1 1]);
   end;
return;

% Included for completeness (usually in own file)
function [fig_pointer_pos, axes_pointer_val] = pointer2d(fig_hndl,axes_hndl)
%
%pointer2d(fig_hndl,axes_hndl)
%
%	Returns the coordinates of the pointer (in pixels)
%	in the desired figure (fig_hndl) and the coordinates
%       in the desired axis (axes coordinates)
%
% Example:
%  figure(1),
%  hold on,
%  for i = 1:1000,
%     [figp,axp]=pointer2d;
%     plot(axp(1),axp(2),'.','EraseMode','none');
%     drawnow;
%  end;
%  hold off

% Rick Hindman - 4/18/01

if (nargin == 0), fig_hndl = gcf; axes_hndl = gca; end;
if (nargin == 1), axes_hndl = get(fig_hndl,'CurrentAxes'); end;

set(fig_hndl,'Units','pixels');

pointer_pos = get(0,'PointerLocation');	%pixels {0,0} lower left
fig_pos = get(fig_hndl,'Position');	%pixels {l,b,w,h}

fig_pointer_pos = pointer_pos - fig_pos([1,2]);
set(fig_hndl,'CurrentPoint',fig_pointer_pos);

if (isempty(axes_hndl)),
	axes_pointer_val = [];
elseif (nargout == 2),
	axes_pointer_line = get(axes_hndl,'CurrentPoint');
	axes_pointer_val = sum(axes_pointer_line)/2;
end;
