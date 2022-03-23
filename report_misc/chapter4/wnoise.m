clear all;clc;close force all;

t.duration = 10;
fs = 48e3;
mu=0;
sd=2/6;

t.vector = 1/fs:1/fs:t.duration;
wn.randn = randn(t.duration*fs , 1);
wn.rand = -1+2*rand(t.duration*fs , 1);
wn.norm = normrnd(mu , sd , t.duration*fs , 1);
wn.randn_norm = normalize(wn.randn, 'range', [-1 .98]);

figure;
sgtitle('Normal distribution (randn)')
subplot(1,2,1)
plot(t.vector , wn.randn)
xlabel('time (s)')
ylabel('Amplitude')
subplot(1,2,2)
histogram(wn.randn,40)
xlabel('Amplitude')
ylabel('Frequency')

set(gcf,'Position',[1   595   977   400])

figure;
sgtitle('Uniform distribution (rand)')
subplot(1,2,1)
plot(t.vector , wn.rand)
xlabel('time (s)')
ylabel('Amplitude')
subplot(1,2,2)
histogram(wn.rand,40)
xlabel('Amplitude')
ylabel('Frequency')

set(gcf,'Position',[978   595   943   400])

figure;
tmsg = sprintf("Normal distribution: mu=%.2f sd=%.2f",mu,sd);
sgtitle(tmsg)
subplot(1,2,1)
plot(t.vector , wn.norm)
xlabel('time (s)')
ylabel('Amplitude')
subplot(1,2,2)
histogram(wn.norm,40)
xlabel('Amplitude')
ylabel('Frequency')

set(gcf,'Position',[1   109   977   400])

figure;
sgtitle({'Normal distribution (randn)','Normalized in a range of [-1 0.98]'})
subplot(1,2,1)
plot(t.vector , wn.randn_norm)
xlabel('time (s)')
ylabel('Amplitude')
subplot(1,2,2)
histogram(wn.randn_norm,40)
xlabel('Amplitude')
ylabel('Frequency')

set(gcf,'Position',[978   109   943   400])

figure;
histogram(wn.randn,40)
hold on
histogram(wn.rand,40)
histogram(wn.norm,40)
histogram(wn.randn_norm,40)