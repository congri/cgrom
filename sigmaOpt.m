function [sgOpt] = sigmaOpt(thetaOpt, nData, nCoarse, sum_XNormSqMean,...
    sumPhiTXmean, sumPhiSq)
%Compute optimal sigma given optimal theta

sigmaSq = (1/(nData*nCoarse))*(sum_XNormSqMean - 2*thetaOpt'*sumPhiTXmean + thetaOpt'*sumPhiSq*thetaOpt);
%avoid imaginary sigma due to numerical precision
if sigmaSq < 0
    sigmaSq = 0;
end
sgOpt = sqrt(sigmaSq);


end

