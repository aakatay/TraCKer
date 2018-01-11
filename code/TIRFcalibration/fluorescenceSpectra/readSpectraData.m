function [fc,fq,fl,ls,sys] = readSpectraData(isBRT,isCy3,dataFN)
% generates MAT data files 
% run in code\TIRFcalibration\fluorescenceSpectra\spectrumData\
    
    % output MAT
    outFN = 'spectrumDataMAT';
    if exist(outFN), rmdir(outFN,'s'); end
    mkdir(outFN);
    PWD = pwd;
    fcPRBmat = [PWD '\spectrumDataMAT\fcPRB'];
    fqPRBmat = [PWD '\spectrumDataMAT\fqPRB'];
    fcBRTmat = [PWD '\spectrumDataMAT\fcBRT'];
    fqBRTmat = [PWD '\spectrumDataMAT\fqBRT'];
    flCYmat  = [PWD '\spectrumDataMAT\flCY'];
    cd(dataFN)

    lsPRB = [488 561 640];
    lsBRT = [473 532 635];
    
    sysBRT.fluoCollect = 0.3; % NA=
    sysPRB.fluoCollect = 0.35; % NA=
    
    %% load Data
    if isBRT
        fcBRT;    %(fc) filter cube 
        fqBRT;    %(fq) quadview/CSU
        ls.X = lsBRT;   %(ls) lasers
        ls.P = [1 1 1]; % laser power
        sys = sysBRT;
    else
        fcPRB;    %(fc) filter cube
        fqPRB;    %(fq) quadview/CSU
        ls.X = lsPRB;   %(ls) lasers 
        ls.P = [1 1 1]; % laser power
    end
    
    if isCy3
        flCY;     %(fl) fluorophores
    else
        flGFP;     %(fl) fluorophores
    end

    %% shape input data to structs
    % fc
    fc.X{1} = fcEx(:,1); % spectra (x-axis)
    fc.X{2} = fcEm(:,1);
    fc.X{3} = fcBs(:,1);
    fc.D{1} = fcEx(:,2); % values (y-axis)
    fc.D{2} = fcEm(:,2);
    fc.D{3} = fcBs(:,2); 

    if isBRT % quadview
        % fqm
        fqm.X{1} = fqm1(:,1); % spectra (x-axis)
        fqm.X{2} = fqm2(:,1);
        fqm.X{3} = fqm3(:,1);
        fqm.D{1} = fqm1(:,2); % values (y-axis)
        fqm.D{2} = fqm2(:,2);
        fqm.D{3} = fqm3(:,2); 
        % fqd
        fqd.X{1} = fqd1(:,1); % spectra (x-axis)
        fqd.X{2} = fqd2(:,1);
        fqd.X{3} = fqd3(:,1);
        fqd.D{1} = fqd1(:,2); % values (y-axis)
        fqd.D{2} = fqd2(:,2);
        fqd.D{3} = fqd3(:,2); 

        fq.fqm = fqm; fq.fqd = fqd;
    else % CSU
        % fqm
        fqm.X{1} = fqm1(:,1); % spectra (x-axis)
        fqm.D{1} = fqm1(:,2); % values (y-axis) 
        fq.fqm = fqm; fq.fqd = [];

    end
    cd(PWD)

    %% KURAL LAB (PRB)
    %filter cube (Chroma 89401 - ET - DAPI/FITC/TRITC/CY5 Quad)
    function fcPRB
        fnX = '89101x.txt';
        fnM = '89101m.txt';
        fnB = '89100bs.txt';
        fq = 'semrockFF01-446-523-600-677-25.txt';
        fcEx = readSpectra(fnX);
        fcEm = readSpectra(fnM);
        fcBs = readSpectra(fnB); 
        fcEx = reshape(fcEx,2,numel(fcEx)/2)';
        fcEm = reshape(fcEm,2,numel(fcEm)/2)';
        fcBs = reshape(fcBs,2,numel(fcBs)/2)';
        save(fcPRBmat,'fcEx','fcEm','fcBs');
    end
    % CSU
    function fqPRB
    %fqm1 = [fcBs(:,1) fcBs(:,2)*0+1]; % none
        fqm1 = readSpectra(fq); % CSU emission filter
        fqm1 = reshape(fqm1,2,numel(fqm1)/2)';

        save(fqPRBmat,'fqm1');
    end

    %% FISHEL LAB (BRT)
    function fcBRT
        fnB = 'semrockDi03-R405-488-532-635-t1-25x36.txt';
        %fnB = 'chroma_zt473-532-633-rpc-uf.txt';
        fcBs = readSpectra(fnB);
        fcBs = reshape(fcBs,2,numel(fcBs)/2)';
        fcEx = [fcBs(:,1) fcBs(:,2)*0+1]; % none
        fcEm = [fcBs(:,1) fcBs(:,2)*0+1]; % none
        save(fcBRTmat,'fcEx','fcEm','fcBs');
    end

    % quadview
    function fqBRT
        fqd1 = readSpectra('QVd1_538semrock.txt');
        fqd1 = reshape(fqd1,2,numel(fqd1)/2)';
        fqd2 = readSpectra('QVd2_640semrock.txt');
        fqd2 = reshape(fqd2,2,numel(fqd2)/2)';
        fqd3 = readSpectra('QVd3_740chroma.txt');
        fqd3 = reshape(fqd3,2,numel(fqd3)/2)';

        fqm1 = readSpectra('QVm1_511-20-25semrock.txt');
        fqm1 = reshape(fqm1,2,numel(fqm1)/2)';
        fqm2 = readSpectra('QVm2_593-46semrock.txt');
        fqm2 = reshape(fqm2,2,numel(fqm2)/2)';
        fqm3 = readSpectra('QVm3_690-50chroma.txt');
        fqm3 = reshape(fqm3,2,numel(fqm3)/2)';
        save(fqBRTmat,'fqd1','fqd2','fqd3','fqm1','fqm2','fqm3');
    end



    %% fluo CY
    function flCY
        cy2fn = 'Alexa Fluor 488.csv';
        cy3fn = 'Cy3.csv';
        cy5fn = 'Cy5.csv';

        cy2 = readSpectra2(cy2fn);
        cy3 = readSpectra2(cy3fn);
        cy5 = readSpectra2(cy5fn);


        % spectra (x-axis)
        fl=[];
        fl.X{1} = cy2(:,1);
        fl.X{2} = cy3(:,1);
        fl.X{3} = cy5(:,1);

        % values (y-axis)
        % excitation
        fl.Dx{1} = cy2(:,2)/100;
        fl.Dx{2} = cy3(:,2)/100;
        fl.Dx{3} = cy5(:,2)/100;
        % emission
        fl.Dm{1} = cy2(:,3)/100;
        fl.Dm{2} = cy3(:,3)/100;
        fl.Dm{3} = cy5(:,3)/100;
        cy3Q = 0.31; % 
        cy3e = 150000; % extinction coeff
        fl.B = [1 1 1]*1e-6; % fluo brightness
        save(flCYmat,'fl');
    end
end
