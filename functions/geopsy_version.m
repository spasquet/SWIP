function version = geopsy_version()
%%% S. Pasquet - V20.04.02

[~,test] = unix('geopsy -app-version');

geopsy = strfind(test,'geopsy ');
version = str2double(test(geopsy(1)+7));