function [yL, yR] = sim_stereo( varargin )

% [1] receiver Structure:
% .type ('0' for receiver0, '1' for receiver with surfaces)
% .struct (microphone and walls structure)
% .loc (receiver location in coordinates [x y z])
% .th (receiver azimuth inclination)
% .phi (receiver z-axis inclination)
% .mic.dmf (microphone distance multiplying factor - dmf)
%
% [2] division/room:
% .size (room size in WxLxH -> [x y z])
% .coeff (room [walls ceiling floor] reflection coefficients)
% .MR (maximum reflections)
% .R (Room() structure)
%
% [3] speaker:
% .loc (speaker location in coordinates [x y z])
%
% [4] inAudio:
% .s (array with all audio samples)
% .fs (audio sampling frequency)
% .fl (lower frequency bound)
% .fh (upper frequency bound)
% 
% ----------------SETTINGS----------------
% [5] hidefig:
% '2' (display figures with virtual reflections)
% '1' (hide figures)
% '0' (display figures)

    if(size(varargin,2) ~= 5)
        error('Wrong number of arguments.');
    end
    % to-do: - comment every line, including the object, etc
    %        - finnish the structure check (numeric, nº elements, etc.)
    
    rec = varargin{1};
    div = varargin{2};
    spk = varargin{3};
    inAudio = varargin{4};
    hidefig = varargin{5};

    rec_ = Receiver(rec.struct, rec.mic.dmf);
    AW = addDivision(div.size, div.coeff);
    
    if rec.type==0
        AM = addReceiver0(rec_, rec.loc, [rec.th rec.phi]);
    elseif rec.type==1
        AM = addReceiver(AW, rec_, rec.loc, [rec.th rec.phi]);
    end

    % show division
    makeFile('headwall', AW, AM, div.MR, inAudio.fs, inAudio.fl, inAudio.fh);
    if hidefig==0
        displayRoom('headwall','HideVS');
    elseif hidefig==2
        displayRoom('headwall');
    end

    % add + draw speaker
    S = addSpk(spk.loc);
    xs=spk.loc(1);
    ys=spk.loc(2);
    zs=spk.loc(3);
    if hidefig==0
        hold on; plot3(xs, ys, zs, 'k*', 'LineWidth', 3); hold off;
    end

    % generate impulse response
    I = impR('headwall', S, div.R); % Compute the impulse response for each microphone

    % generate L + R channel (2 channels)
    yL = fftfilt(I(1,:),inAudio.s);
    yR = fftfilt(I(2,:),inAudio.s);
end

