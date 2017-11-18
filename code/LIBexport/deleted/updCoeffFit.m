% reads the parameters from the sliders and sets the coefficent values
% accordingly
                for i = 1:numCoeffFitParam
                    coeffParam(i) = get(hImgCoeffParam(i),'Value');
                    set(hImgCoeffParamText(i),'String',coeffParam(i)*Coeff);
                end
                xPoly = [1 [1:numCoeffFitParam]*coeffFitParamCoverage];
                xPoly(end) = Frames-frstFrm2+1;
                CoeffFit = linearfit(xPoly,[1 coeffParam]);
                set(coeffPlot,'YData',CoeffFit);