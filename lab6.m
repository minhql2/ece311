%% Ideal Filters

%% lowpass 20
N = 20;
n = linspace(-N,N,2*N);
wc = pi/2;
h1 = wc / pi * sinc(wc / pi * n);
H1 = fft(h1);
H1 = fftshift(H1);
N = 2*N;
w = fftshift((0:N-1)/N * 2 * pi);
w(1:N/2) = w(1:N/2) - 2*pi;
figure(1);
subplot(2,1,1);
plot(n, h1);
subplot(2,1,2);
stem(w,abs(H1));

%% highpass 20
N = 20;
y = dirac(n);
y(21) = 1;
wc = pi/2;
h2 = wc / pi * sinc(wc / pi * n);
H2 = fft(h2);
Y = fft(y);
H2 = Y - H2;
H2 = fftshift(H2);
figure(2);
subplot(2,1,1);
plot(n, h2);
subplot(2,1,2);
stem(w,abs(H2));

%% bandpass 20
N=20;
wo = pi/2;
wc = pi/4;
h3 = cos(wo * n).* wc / pi .* sinc(wc / pi * n);
H3 = fft(h3);
H3 = fftshift(H3);
figure(3);
subplot(2,1,1);
plot(n, h3);
subplot(2,1,2);
stem(w,abs(H3));

%% lowpass 40
N = 40;
n = linspace(-N,N,2*N);
wc = pi/2;
h1 = wc / pi * sinc(wc / pi * n);
H1 = fft(h1);
H1 = fftshift(H1);
N = 2*N;
w = fftshift((0:N-1)/N * 2 * pi);
w(1:N/2) = w(1:N/2) - 2*pi;
figure(1);
subplot(2,1,1);
plot(n, h1);
subplot(2,1,2);
stem(w,abs(H1));

%% highpass 40
y = dirac(n);
y(41) = 1;
wc = pi/2;
h2 = wc / pi * sinc(wc / pi * n);
H2 = fft(h2);
Y = fft(y);
H2 = Y - H2;
H2 = fftshift(H2);
figure(2);
subplot(2,1,1);
plot(n, h2);
subplot(2,1,2);
stem(w,abs(H2));

%% bandpass 40
wo = pi/2;
wc = pi/4;
h3 = cos(wo * n).* wc / pi .* sinc(wc / pi * n);
H3 = fft(h3);
H3 = fftshift(H3);
figure(3);
subplot(2,1,1);
plot(n, h3);
subplot(2,1,2);
stem(w,abs(H3));



%% Report Item 2
load impulseresponse.mat;
n = linspace(0,73,74);
N = 512;
w = fftshift((0:N-1)/N * 2 * pi);
w(1:N/2) = w(1:N/2) - 2*pi;
H = fft(h, 512);
H = fftshift(H);
figure(4);
plot(w, mag2db(abs(H)));


%% Report Item 3
N = 25;
w = fftshift((0:N-1)/N * 2 * pi);
w(1:12) = w(1:12) - 2*pi;
n = linspace(0,24,25);
wc = pi / 3;
M= (N-1)/2;
G = rectangularPulse(-wc,wc,w) .* exp(-1j * M * w);
G = ifftshift(G);
g = ifft(G);

window = hamming(25);
window = transpose(window);
h = window .* g;
H = fft(h,512);
H = fftshift(H);
N=512;
w = fftshift((0:N-1)/N * 2 * pi);
w(1:N/2) = w(1:N/2) - 2*pi;
figure(6);
subplot(2,1,1);
plot(n,h);
subplot(2,1,2);
plot(w, mag2db(abs(H)));
