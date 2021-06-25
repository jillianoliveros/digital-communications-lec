%Task 1

% note: the imported text file can be the short version 
% or the longer version. kindly change the file name below.
% text_file_long_ver.txt (long version)
% text_file_short_ver.txt (short version)

clc
clear
s=importdata('text_file_short_ver.txt');
str=[s(:)];
str=[str{:}]; %walang /n sa appendix kaya nag eerror.
concatwSC=['ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz,. ;']; %appendix  
matching=sum(ismember(str,concatwSC)); %
i=length(concatwSC)+1;%init ng while statement
x=1;

while i > 1
   i=i-1;
  freq(i,1)=sum(ismember(str,concatwSC(i)))/matching; % pang solve ng frequency ng kada letters
end
concatwSC=num2cell(concatwSC);%convert kada cell kasi double yung data type ng concatwSC
Appendix=[transpose(concatwSC) num2cell(freq)]; %para makita yung frequency ng kada symbol
[huff,man]=huffmandict(concatwSC,freq);%para sa tree diagram 
enco=huffmanenco(str,huff);%pang encode
%mess=huffmandeco(enco,huff); %ilalagay na yung mga nadecode na latters
                                %kaso naka cell lang sya
%message=cell2mat(mess);    
 %--------Entropy-------%

for x = 1 : length(freq)

 Entrop(x) = freq(x) * log(1/freq(x))/log(2) ;
end
Entropy=nansum(Entrop);
fprintf('The Entropy is %d bits/symbol\n', Entropy)
fprintf('The Average Length is %d bits/symbol\n', man)


%Task 2

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

%Oversampling
fs = length(x);
s = 10;
oversampled = repelem(x,s);
figure;
stem(oversampled);

%Pulse-shaping filter 
Beta = 0.95; %pick the parameters
span = 1; % filter span in symbols
sps = s; %ovesampling factor
rcs = rcosdesign(Beta,span,sps,'normal');
upf = upfirdn(oversampled,rcs,span);
t = (0:(length(upf)-1))/fs;
figure;
plot(t,upf);


%Task 3

in = input('Choose a modulation scheme (0 for BPSK, 1 for BASK): '); 
 
% map the continuous baseband waveform from task 2 into 
% bits suitable for transmission
binary = [];
if (in == 0)
    % bpsk implementation
    for (i = 1:1:length(upf))
        if (upf(i)>0)
            bi = 1;
        elseif (upf(i)<0)
            bi = -1;
        end
        binary = [binary bi];
    end
elseif (in == 1)
    % bask implementation
    for (i = 1:1:length(upf))
        if (upf(i)>0)
            bi = 1;
        elseif (upf(i)<0)
            bi = 0;
        end
        binary = [binary bi];
    end
end
 
subplot(3,1,1)
plot(t,binary)
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('Binary Information');
 
% bpsk/bask modulation
fc = 3; % frequency of the carrier
carrier = sin(2*pi*fc*t); % equation of the carrier signal
mod = binary.*carrier; % generate the modulated signal by multiplying it with the carrier
subplot(3,1,2)
plot(t,carrier); % carrier signal
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('Carrier signal');
subplot(3,1,3)
plot(t,mod); % bpsk/bask modulated signal
xlabel('Time (seconds)');
ylabel('Amplitude (volts)');
title('Modulated signal');
 
% signal constellation plot
figure
set(gcf,'color','w')
plot(real(binary),imag(binary),'o') % constellation diagram
grid on; box on;
title('Constellation Diagram')
xlim([-1.5 1.5]); ylim([-1.5 1.5]);


%Task 4

% Energy-per-bit-to-noise ratio (range)
EbN0_dB = -10:2:10; % Eb/N0 in dB
EbN0 = 10.^(EbN0_dB/10); % convert Eb/No to linear scale

Eb = sps * norm(mod)^2/length(mod); % energy
N0 =  Eb./EbN0; % noise

for i = 1:length(EbN0)
    
    mod_2 = mod;
    if (iscolumn(mod))
        mod = mod.';
    end

    mod_l = size(mod);
    if (isreal(mod))
        n_var = (N0(i))./2; % noise variance
        s_dev = sqrt(n_var); % standard deviation for AWGN Noise when 'mod' is real
        n = s_dev * randn(mod_l); % computed noise
    else
        n_var = (N0(i))./2; % noise variance
        s_dev = sqrt(n_var/2); % standard deviation for AWGN Noise when 'mod' is complex
        n = s_dev * (randn(mod_l) + 1i.*randn(mod_l)); % computed noise
    end
    
    if (iscolumn(mod_2))
        n = n.';
    end

    r = mod + n;

figure
plot(n(1:1000))
grid on;
title ('Noise')
figure
set(gcf,'color','w')
plot(real(r),imag(r),'o')
grid on; box on;
title('Constellation Diagram')
xlim([-1.5 1.5]); ylim([-1.5 1.5]);


% Task 5

    %Front-End Receiver
    cs = sin(2*pi*fc*t); %carrier signal
    rx = r(i).*cs(1:length(r)); %signal mixing

    %Integrator
    lb = length(bits);
    rx_integ = conv(rx,ones(1,lb)) * 2/lb; 
    % Sample implementation of integrator

    %Downsampling
    rx_down = downsample(rx_integ,lb); %bult-in

    %Thresholding
    for k=1:length(rx_down)
        if rx_down(k)<=0 
          b_hat(k)=0;
        else
          b_hat(k)=1;

        end
    end

    %BER computation (experimental)
    E = abs(bits-b_hat(2:length(bits)+1)); 
    eBER(i) = sum(E)/length(bits);

end

%BER computation (theoretical)
tBER = (1/2)*erfc(EbN0/sqrt(2));

%Plot
figure
semilogy(EbN0_dB,tBER,'r-')
hold on
semilogy (EbN0_dB,eBER,'go')
grid on
legend('Theoretical BER', 'Experimental BER');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve')
