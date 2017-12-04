load('rtTraCKerTraceTIME')
figure;

    %% rtTraCKerTrace
    subplot(1,2,1)
    area([tprm(:,1)';tprm(:,2)';tprm(:,3)';tprm(:,4)';tprm(:,5)';tprm(:,6)']')
    hold on
    TLOOP = tloop(2:end)-tloop(1:end-1);
    TLOOP(end) = TLOOP(end-1);
    plot(TLOOP)
    hold off
    title(sprintf('mean time : %.03f',mean(tprm(:,1))));
    legend('factorfast','trace','clearIMG','genImg','intro','save','total')
    ylim([0 1])
    
    
    %% factorFast
    subplot(1,2,2)
    plot(ffT);
    legend('factorFastTotal','isMember')
    title(sprintf('mean time : %.03f',mean(ffT(:,2))));
    ylim([0 1])
    