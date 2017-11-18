function [TIRFangle] = calcTIRFangle(TIRFstagePosEPI,TIRFstagePosTIRF)
    coeffTIRF_EPI_conversion = 0;
    TIRFangle = (TIRFstagePosEPI-TIRFstagePosTIRF)*coeffTIRF_EPI_conversion;
end