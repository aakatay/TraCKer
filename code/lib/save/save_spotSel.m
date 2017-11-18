save('spotSel','tP2G_1','tP2G_2','tol','tol0',...
'sg0','sg0min','sg0max','sg','sr','intP','minmaxIntTol');

cfg.fit.gaus = struct;    
cfg.fit.gaus.tP2G_1 = tP2G_1;
cfg.fit.gaus.tP2G_2 = tP2G_2;
cfg.fit.gaus.tol = tol;
cfg.fit.gaus.tol0 = tol0;

cfg.fit.gaus.sg0 = sg0;
cfg.fit.gaus.sg0max = sg0max;
cfg.fit.gaus.sg0min = sg0min;
cfg.fit.gaus.sg = sg;
cfg.fit.gaus.sr = sr;
cfg.fit.gaus.intP = intP;