close all;
clear;
clc;
%Manchester Line Encoding.
bits = enco;
bitrate = 1;
n = 1000;
T = length(bits)/bitrate;
N = n*length(bits);
dt = T/N;
t = 0:dt:T;
x = zeros(1,length(t));
for i=1:length(bits)
  if bits(i)==1
    x((i-1)*n+1:(i-1)*n+n/2) = 1;
    x((i-1)*n+n/2:i*n) = -1;
  else
    x((i-1)*n+1:(i-1)*n+n/2) = -1;
    x((i-1)*n+n/2:i*n) = 1;
  end
end
 plot(t, x, 'Linewidth', 3);

% counter = 0;
% for i = 1:length(t)
%   if t(i)>counter
%     counter = counter + 1;
%     if x(i)>0
%       result(counter) = x(i);
%     else result(counter) = 0;
%     end
%   end
% end
% disp('Manchester Decoding:');
% disp(result);

%Oversampling
fs = length(x);
% tfs = 0:1/fs:1-1/fs;
s = 10;
oversampled = repelem(x,s);
% nfs = 0:0.1:length(x);  %new sampling rate
figure;
stem(oversampled);
%pulse shaping filter 
Beta = 0.95; %pick the parameters
span = 1; % filter span in symbols
sps = s; %ovesampling factor
rcs = rcosdesign(Beta,span,sps,'normal');
upf = upfirdn(rcs,oversampled,span);
t = (0:(length(upf)-1))/fs;
figure;
plot(t,upf);

binary = [];
for (i = 1:1:length(upf))
    if (upf(i)>0)
        bi = 1;
    else
        bi = -1;
    end
    binary = [binary bi];
end
figure;
plot(t,binary)
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('Binary Information');

%bpsk modulation
f = 3; % frequency of the carrier
carrier = sin(2*pi*f*t);
mod = binary.*carrier; % generate the modulated signal by multiplying it with the carrier
figure;
plot(carrier); % carrier signal
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('Carrier signal');
figure;
plot(mod); % bpsk modulated signal
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('BPSK Modulated signal');

% signal constellation plot
figure
set(gcf,'color','w')
plot(real(binary),imag(binary),'o')
grid on; box on;
title('Constellation Diagram')
xlim([-1.5 1.5]); ylim([-1.5 1.5]);