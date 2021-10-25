function [AOA,DELAY_t] = detect_az2(CH, CR, C, D_X)

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
    [MAX, IDX]=         max(CH_CORR);
    % INTERPOLATION
    if IDX > 1 && IDX < length(CH_CORR)
        OFFSET_i=       (CH_CORR(IDX-1) - CH_CORR(IDX+1))/...
                        (2*CH_CORR(IDX-1) + 2*CH_CORR(IDX+1) - 4*MAX);
        IDX = IDX + OFFSET_i;
    end
    % END OF INTERPOLATION
    DELAY_i=            (IDX-1) - CR;
    DELAY_t=            DELAY_i/CH.fs;
    arg_=               C*DELAY_t/(D_X);
    arg_(arg_>1)=       1;
    arg_(arg_<-1)=      -1;
    AOA=                acosd(arg_);
    
end

