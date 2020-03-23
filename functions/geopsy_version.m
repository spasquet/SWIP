function version = geopsy_version()

[~,test] = unix('dinver -version');

geopsy = strfind(test,'geopsypack-');
version = str2double(test(geopsy(1)+11));