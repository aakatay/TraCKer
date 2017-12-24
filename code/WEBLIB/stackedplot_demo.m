clear; close all; drawnow

N = 4;
mu = 2*rand(1,N);
sigma = [1 0.1 10 1];%rand(1,N);

t = 1:100;
x = bsxfun(@plus,randn(100,4)*diag(sigma),mu);

subplot(221)
plot(t,x)
title('plot(t,x)')
subplot(222)
stackedplot(t,x,'PlotSpans','dynamic','YAxisLocation','right')
title('stackedplot(t,x,''PlotSpans'',''dynamic'',''YAxisLocation'',''right'')')
subplot(223)
stackedplot(t,x,'PlotSpans','equal','DataTickMode','baseline','BaseValue',0)
title(sprintf('stackedplot(t,x,''PlotSpans'',''equal'',...\n''DataTickMode'',''baseline'',''BaseValue'',0)'))
subplot(224)
stackedplot(x,t,'PlotScaling',[1 10 0.1 1],'DataTickMode','all','PlotSpacing',0.5)
title(sprintf('stackedplot(x,t,''PlotScaling'',[1 10 0.1 1],...\n''DataTickMode'',''all'',''PlotSpacing'',0.5)'))
set(gca,'FontSize',8)

t = 1:100;
x = bsxfun(@plus,rand(100,4)*diag(sigma),mu);

figure
subplot(121)
semilogy(t,x)
title('semilogy(t,x)')
subplot(122)
stackedplot(t,x,'YAxisScale','log','LabelLocation','baseline','ShowBaseLine',false)
title('stackedplot(t,x,''YAxisScale'',''log'',''LabelLocation'',''baseline'',''ShowBaseLine'',false)')
