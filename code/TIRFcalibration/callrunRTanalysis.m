try runRTanalysis
catch ME
    deleteALLparallelPools
    gcp
    fclose('all')
    cd('E:\MATLAB\TIRFcalibration\data\Ata01_14_realtime_interface')
    rethrow ME
end