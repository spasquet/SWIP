function cmap = color_2_white(col,ncol)

col1 = linspace(col(1),1,ncol);
col2 = linspace(col(2),1,ncol);
col3 = linspace(col(3),1,ncol);

cmap = [col1' col2' col3'];
