function swip_version
%%% S. Pasquet - V16.5.4
% Give SWIP current version

Ascript=which('A_SWIPdisp_script');
Astruct=dir(Ascript);
if isempty(Astruct)==0
    Adate=Astruct.date;
    fprintf(['\n A_SWIPdisp_script installed in ',fileparts(Ascript)]);
    fprintf(['\n Last modification: ',Adate,'\n']);
else
    fprintf('\n A_SWIPdisp_script not found\n');
end

Bscript=which('B_SWIPparam_script');
Bstruct=dir(Bscript);
if isempty(Bstruct)==0
    Bdate=Bstruct.date;
    fprintf(['\n B_SWIPparam_script installed in ',fileparts(Bscript)]);
    fprintf(['\n Last modification: ',Bdate,'\n']);
else
    fprintf('\n B_SWIPparam_script not found\n');
end

Cscript=which('C_SWIPinv_script');
Cstruct=dir(Cscript);
if isempty(Cstruct)==0
    Cdate=Cstruct.date;
    fprintf(['\n C_SWIPinv_script installed in ',fileparts(Cscript)]);
    fprintf(['\n Last modification: ',Cdate,'\n']);
else
    fprintf('\n C_SWIPinv_script not found\n');
end

D1script=which('D1_SWIPmod1d_script');
D1struct=dir(D1script);
if isempty(D1struct)==0
    D1date=D1struct.date;
    fprintf(['\n D1_SWIPmod1d_script installed in ',fileparts(D1script)]);
    fprintf(['\n Last modification: ',D1date,'\n']);
else
    fprintf('\n D1_SWIPmod1d_script not found\n');
end

D2script=which('D2_SWIPmod2d_script');
D2struct=dir(D2script);
if isempty(D2struct)==0
    D2date=D2struct.date;
    fprintf(['\n D2_SWIPmod2d_script installed in ',fileparts(D2script)]);
    fprintf(['\n Last modification: ',D2date,'\n']);
else
    fprintf('\n D2_SWIPmod2d_script not found\n');
end

defscript=which('SWIP_defaultsettings');
defstruct=dir(defscript);
if isempty(defstruct)==0
    defdate=defstruct.date;
    fprintf(['\n SWIP_defaultsettings installed in ',fileparts(defscript)]);
    fprintf(['\n Last modification: ',defdate,'\n']);
else
    fprintf('\n SWIP_defaultsettings not found\n');
end