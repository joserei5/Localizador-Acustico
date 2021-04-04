classdef Room < handle
    % Ultrasonic Propagation Object to be used in the impR function.
    
    properties
        T = 22;     % Room Temperature
        P = 1;      % Room Pressure
        H = 50;     % Room Humidity
        fs = 160e3; % Sampling frequency
        fl = [];    % Lower bound
        fh = [];    % Upper Bound for frequency optimization
        
        L_d = 32;   % Delay filter order
        M   = 100;  % Number of iterations
        L_a = 32;   % Attenuation filter order
        I   = 10;   % FFT iterpolation factor (for sound attenuation)
        
        dref = 1;   % reference distance        
    end
    
    properties(SetAccess = private, GetAccess = private)
        c = [];
        a = [];
        f = [];
    end
    
    methods
        % Class Constructor
        function obj = USC_TC1(fs)
            if(nargin==1)
                obj.fs = fs;
            else
                refresh(obj);
            end
        end
        
        function obj = set.fs(obj,value)
            obj.fs = value;
            refresh(obj);
        end
        
        function obj = set.T(obj,value)
            obj.T = value;
            refresh(obj);
        end
        
        function obj = set.P(obj,value)
            obj.P = value;
            refresh(obj);
        end
        
        function obj = set.H(obj,value)
            obj.H = value;
            refresh(obj);
        end
        
        
        function [h,p] = impulse(obj,d)
            
            % Fractional delay
            
            dn = d*obj.fs/obj.c;
            
            [hd,pd] = Delay(obj,dn);
            
            % Attenuation with frequency
            
            [ha,pa] = DistanceAtt(obj,d);
            
            h = conv(hd,ha);
            p = pd+pa;
            
        end
        
    end
    
    methods (Access = private)
        
        function obj = refresh(obj)
            
            obj.f = 0:2*obj.fs/(obj.L_a*obj.I):obj.fs/2; % half frequency vector
            
            obj.a = SoundAtt(obj.P, obj.T, obj.H, obj.f);   %   absorption attenuation for each frequency
            
            obj.c = sqrt(1.098395978344934e+05*(1+obj.T/273));
        end
        
        function [h,p] = Delay(obj,dn)
            
            L = obj.L_d+1;            % Fractional delay filter coefficients
            L2 = floor(L/2);              % Integer delay of the filter
            
            di = floor(dn);
            df = dn-di;
            
            p = di-L2;
            
            h = ones(1,L);                %filter coefficients
            D=L2+df;
            
            %compute the filter coefficients
            for n=0:obj.L_d
                for k=0:obj.L_d
                    if (k~=n)
                        h(n+1) = h(n+1)*(D-k)/(n-k);
                    end  % if
                end % for k
            end % for n
            
        end
        
        function [h,p] = DistanceAtt(obj,d)
            
            at = obj.dref*10.^(-obj.a.*(d-obj.dref)./20)/d;
            
            S = [at at(end-1:-1:2)]; %attenuation spectrum
            
            N = length(S);
            
            D = floor(N/2);
            
            pa = -D;
            
            n = 0:N-1;
            
            SD = S.* exp(-1i*2*pi*n*D./N);
            
            sd = ifft(SD);
            
            sd = real(sd); % remove some aproximation error
            
            L = obj.L_a+1;
            
            px = floor((N-L)/2);
            
            mask = ([zeros(1,px) ones(1,L) zeros(1,ceil((N-L)/2))]);
            
            sd = sd.*mask;
            
            %%% with interative method
            if(~isempty(obj.fl) && ~isempty(obj.fh))
                N2 = length(at);
                Nci = floor((N2-1)*2*obj.fl/obj.fs);
                Nce = ceil((N2-1)*2*obj.fh/obj.fs+1);
                
                % freq. window
                nm = logical([zeros(1,Nci) ones(1,Nce-Nci) zeros(1,N-2*Nce+1) ones(1,Nce-Nci) zeros(1,Nci-1)]);
                if(Nci==0)
                    nm(end)=[];
                end
                
                
                for m=1:obj.M
                    SDi = fft(sd);
                    SDi(nm) = SD(nm);
                    sd = ifft(SDi);
                    sd = sd.*mask;
                    sd = real(sd); % remove some aproximation error
                end
                
            end
            
            h = sd(floor((N-L)/2)+1:floor((N-L)/2)+L);
            
            p = pa+px;
            
        end
        
        
    end
    
end