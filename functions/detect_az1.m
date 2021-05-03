function AOA = detect_az1(CH, CR, C, D_X)

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

    CH_CORR=            xcorr(CH.L,CH.R,CR);
    [~, CORR_i]=        max(abs(CH_CORR));
    DELAY_i=            (CORR_i-1) - CR;
    DELAY_t=            DELAY_i/CH.fs;
    arg_=               C*DELAY_t/(D_X);
    arg_(arg_>1)=       1;
    arg_(arg_<-1)=      -1;
    AOA=                acosd(arg_);
    
end

