function [sectDataX,sectDataY,sectFitX,sectFitY] = getCrossSect(IdatSHOW,Ifit2SHOW_2,posSect,intPeak)
    
    p1X = round(posSect(1));
    p1Y = round(posSect(2));
    p2X = round(posSect(3));
    p2Y = round(posSect(4));
    sectDataX = IdatSHOW(:,p1X);
    sectDataY = IdatSHOW(p1Y,:);
    sectFitX = Ifit2SHOW_2(:,p2X)*intPeak;
    sectFitY = Ifit2SHOW_2(p2Y,:)*intPeak;
end