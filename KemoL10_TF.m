function [h, p] = KemoL10_TF(theta,phi,fs,fi,fe)

load Kemo_L10
load Kemo_L10_RAD

df = 50;    % Frequency resolution
K = 260e-6; % Original delay
K = round(K*fs);  
p = 20;     % Impulse response artificial delay
De = K-p;   % Final delay correction

L = 256;  % Final impulse length
Mi = 10;  % Number interactions

Ni = fi/df+1;
Ne = fe/df+1;

f = 0:df:fs-df;

N = length(f);

D = floor(N/2);
n = 0:N-1;
E = exp(-1i*2*pi*n*D./N);

% freq. window
nm = logical([zeros(1,Ni-1) ones(1,Ne-Ni+1) zeros(1,N-2*Ne+1) ones(1,Ne-Ni+1) zeros(1,Ni-2)]);
% impulse mask
mask = ([zeros(1,D+De) ones(1,L) zeros(1,N-D-L-De)]);

[x,y,z] = sph2cart(theta,phi,1);
a = acos(x/sqrt(x.^2+y.^2+z.^2)); % get the angle

Xri = interp2(L10R.a,L10R.f,L10R.X,a,f(Ni:Ne)','spline').';
Xsi = interp1(L10.f,L10.p,f(Ni:Ne),'spline');

Xi = Xri.*Xsi;

Xii = [zeros(1,Ni-1) Xi zeros(1,N-2*Ne+1) conj(Xi(end:-1:1)) zeros(1,Ni-2)];
%keyboard
Xii = Xii.*E;

x = ifft(Xii);
for m=1:Mi
    X = fft(x);
    X(nm) = Xii(nm);
    x = ifft(X);
    x = x.*mask;
    x = real(x); % remove some aproximation error
end

h = x(logical(mask));
