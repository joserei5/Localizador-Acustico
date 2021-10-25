function [AOA,DELAY_t] = detect_az3(CH, CR, C, D_X)

% [1] CH:
% .L (left channel)
% .R (right channel)
% .fs (sampling frequency)
%
% [2] CR:
% correlation range
%
% [3] C:
% sound velocity (m/s)
%
% [4] D_X:
% receiver's microphone spacing distance (m)
%
% [OUTPUT]
% AOA = experimental angle-of-arrival
% DELAY_t = experimental delay (+/- ms)

    CH_CORR=            abs(xcorr(CH.L,CH.R,CR));
    CR_RANGE=           -CR:0.1:CR;
    CH_CORR=            interp1(-CR:CR,CH_CORR,CR_RANGE,'spline');
%     CH_CORR=            spline(-CR:CR,CH_CORR,CR_RANGE);
    [MAX, IDX]=         max(CH_CORR);
    DELAY_i = CR_RANGE(IDX);
    DELAY_t=            DELAY_i/CH.fs;
    arg_=               C*DELAY_t/(D_X);
    arg_(arg_>1)=       1;
    arg_(arg_<-1)=      -1;
    AOA=                acosd(arg_);
    
end

