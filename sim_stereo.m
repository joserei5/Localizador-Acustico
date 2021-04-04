function [yL, yR] = sim_stereo( varargin )
    %{
        [1] structure = object to simulate
        [2] division = division coord.      = [x y z]
        [3] reclocation = receiver coord.   = [x y z]
        [4] recdir = receiver directions    = [theta phi]
        [5] spklocation = speaker location  = [x y z]
        [6] y  = mono-channel audio sample
        [7] fs = sampling frequency
        ------------------------------------------------
        SETTINGS:
        [8] hidefig = hide figures (1) / display figures (0)
        ------------------------------------------------
        OPTIONAL:
        [9] dmf = structure distance multiplying factor
    %}

    if(size(varargin,2) < 8)
        error('Wrong number of arguments.');
    end
    % to-do: - comment every line, including the object, etc
    %        - finnish the structure check (numeric, nº elements, etc.)
    
    structure = varargin{1};
    division = varargin{2};
    reclocation = varargin{3};
    recdir = varargin{4};
    spklocation = varargin{5};
    y = varargin{6};
    fs = varargin{7};
    hidefig = varargin{8};
    
    if size(varargin,2) > 8
        dmf = varargin{9};
    else
        dmf = 1;
    end

    rec = Receiver(structure, dmf);
    AW = addDivision(division,[0 0 0]);
    AM = addReceiver0(rec, reclocation, recdir);

    MR = 0;     % max reflections (order)
    fl = 100;   % lower frequency bound
    fh = 20e3;  % upper frequency bound (<=80 kHz due to @KemoL10_TF)

    % show division
    makeFile('headwall', AW, AM, MR, fs, fl, fh);
    if hidefig==0
        displayRoom('headwall','HideVS');
    end

    % add + draw speaker
    S = addSpk(spklocation);
    xs=spklocation(1);
    ys=spklocation(2);
    zs=spklocation(3);
    if hidefig==0
        hold on; plot3(xs, ys, zs, 'k*', 'LineWidth', 3); hold off;
    end

    % generate room + impulse response
    R = Room();
    R.T = 25;  % temperatura ºC
    R.H = 50;  % humidade %
    R.P = 1.01;% pressure atm

    I = impR('headwall', S, R); % Compute the impulse response for each microphone
%     figure
%     plot(I')
    % generate L + R channel (2 channels)
    yL = fftfilt(I(1,:),y);
    yR = fftfilt(I(2,:),y);
end

