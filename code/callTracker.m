clear; isCallTracker = 1;
fn = dir('*.tif');

for i=1:size(fn,1)
    fname = fn(i).name;
    display(sprintf('==>%s\n',fname));
    save('i','i','fn');
    save('fname','fname')
    TraCKer_3D_w_ZcolorPlot_deep
    load i;
end