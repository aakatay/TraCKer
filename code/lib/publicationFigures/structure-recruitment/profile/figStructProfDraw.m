% draws profile for profile selection (in selectStruct)

%% profiles

[~,~,profInt] = improfile(repelem(A,4,4),[P1(1),P1(2)],[P1(3),P1(4)]);
[~,~,profRec] = improfile(R,[P1(1),P1(2)],[P1(3),P1(4)]);


figure(101) % profile
sz0 = 1.5; sz1 = 1.5; 
nm = numel(profRec);
v=1:nm;
%plot( profInt/max(profInt)*max(profRec) ,'k','LineWidth',sz0);

hp=plotyy(v,profRec,v,profInt); 
hpc(1)=hp(1).Children;
hpc(2)=hp(2).Children;
set(hpc(1),'LineWidth',sz0);
set(hpc(2),'LineWidth',sz0);
set(hp(1),'LineWidth',sz0);
set(hp(2),'LineWidth',sz0);

%legend('Intensity', 'Recruitments')
title('profile plots')

ylabel(hp(1),'number of recruitments')
ylabel(hp(2),'intensity')
grid minor;
set(gcf,'color','w');