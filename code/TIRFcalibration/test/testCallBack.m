function testCallBack
    close all
    F = findall(0, 'HandleVisibility', 'off'); delete(F)
    fig= uifigure;
    x1 = 50;
    x2 = 90;
    y1 = [1:5]*50+40;
    btn = uibutton(fig,'Position',[20 30 100 20],'Text','Start');

    lmp1 = uilamp(fig,'Position',[x1 y1(1) 20 20],'Color',[1 0 0]);
    lmp2 = uilamp(fig,'Position',[x1 y1(2) 20 20],'Color',[1 0 0]);
    lmp3 = uilamp(fig,'Position',[x1 y1(3) 20 20],'Color',[1 0 0]);
    lmp4 = uilamp(fig,'Position',[x1 y1(4) 20 20],'Color',[1 0 0]);
    lmp5 = uilamp(fig,'Position',[x1 y1(5) 20 20],'Color',[1 0 0]);


    label1 = uilabel(fig,'Position',[x2 y1(1) 100 15],'Text','WAmean');
    label2 = uilabel(fig,'Position',[x2 y1(2) 100 15],'Text','detectThresh');
    label3 = uilabel(fig,'Position',[x2 y1(3) 100 15],'Text','TraCKerPos');
    label4 = uilabel(fig,'Position',[x2 y1(4) 100 15],'Text','TraCKerTrace');
    label5 = uilabel(fig,'Position',[x2 y1(5) 100 15],'Text','trackSNR');

    %% verify running codes
    set(lmp1,'Color',[0 1 0]);
    callBackFunc1(lmp3)
    %%
    %WAmean;
    p = gcp();
    f=parfeval(p,@callBackFunc1,0,lmp2);


end
