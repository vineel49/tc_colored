% Optimum Predictive Iterative Decoding of Turbo Codes in Coloured Gaussian Noise

close all
clear all
clc
%---------------- SIMULATION PARAMETERS ------------------------------------
SNR_dB = -10; % SNR per bit in dB (in logarithmic scale)
sim_runs = 1*(10^1); % simulation runs
frame_size = 1024; % frame size
num_bit = 0.5*frame_size; % number of data bits (overall rate is 1/2)
SNR = 10^(0.1*SNR_dB); % SNR per bit in linear scale
noise_var_1D = 2*2/(2*SNR); % 1D noise variance
%--------------------------------------------------------------------------
%    Generator polynomial of the component encoders
gen_poly = ldiv2([1 0 1],[1 1 1],num_bit); % using long division method

%  Interleaver and deinterleaver mapping of the turbo code 
intr_map = randperm(num_bit);
deintr_map = deintrlv((1:num_bit),intr_map);
%--------------------------------------------------------------------------
% IIR FILTER PARAMETERS (USED TO GENERATE COLOURED NOISE)
a = 0.99;
B = sqrt(1-a^2); % INPUT COEFFICIENTS IN THE DIFFERENCE EQUATION
A = [1 -a]; % OUTPUT COEFFICIENTS IN THE DIFFERENCE EQUATION
AUTOCORR_SEQ = [noise_var_1D; a*noise_var_1D]; % AUTOCORRELATION OF COLOURED NOISE

% GENERATE PREDICTION FILTER COEFFICIENTS USING LEVINSON DURBIN ALGORITHM
[pred_coef,pred_var] = Gen_Coef(AUTOCORR_SEQ,1);

%--------------------------------------------------------------------------
C_Ber = 0; % channel erros
tic()
%--------------------------------------------------------------------------
for frame_cnt = 1:sim_runs
%                           TRANSMITTER
%Source
a = randi([0 1],1,num_bit); % data

% Turbo encoder
% component encoder 1
b1 = zeros(1,2*num_bit); % encoder 1 output initialization
b1(1:2:end) = a; % systematic bit
temp1 = mod(conv(gen_poly,a),2); % linear convolution with the generator polynomial
b1(2:2:end) = temp1(1:num_bit); % parity bit
% component encoder 2
b2 = zeros(1,2*num_bit); % encoder 2 output initialization
b2(1:2:end) = a(intr_map); % systematic bit
temp2 = mod(conv(gen_poly,b2(1:2:end)),2); % linear convolution with the generator polynomial
b2(2:2:end) = temp2(1:num_bit); % parity bit

% QPSK mapping (according to the set partitioning principles)
mod_sig1 = 1-2*b1(1:2:end) + 1i*(1-2*b1(2:2:end));
mod_sig2 = 1-2*b2(1:2:end) + 1i*(1-2*b2(2:2:end));
mod_sig = [mod_sig1 mod_sig2];

%--------------------------------------------------------------------------
%                            CHANNEL   
% AWGN
white_noise = sqrt(noise_var_1D)*randn(1,10000+frame_size)+1i*sqrt(noise_var_1D)*randn(1,10000+frame_size); 

colored_noise = filter(B,A,white_noise);
colored_noise(1:10000) = []; % DISCARDING TRANSIENT SAMPLES
Chan_Op = mod_sig + colored_noise(1:frame_size); % Chan_Op stands for channel output
%--------------------------------------------------------------------------
%                          RECEIVER 

% Branch metrices for the BCJR
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
Dist = zeros(16,frame_size);
Dist(1,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(1))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(1)))).^2;
Dist(2,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(4))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(1)))).^2;

Dist(3,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(1))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(4)))).^2;
Dist(4,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(4))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(4)))).^2;

Dist(5,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(1))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(2)))).^2;
Dist(6,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(4))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(2)))).^2;

Dist(7,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(1))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(3)))).^2;
Dist(8,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(4))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(3)))).^2;

Dist(9,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(2))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(1)))).^2;
Dist(10,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(3))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(1)))).^2;

Dist(11,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(2))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(4)))).^2;
Dist(12,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(3))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(4)))).^2;

Dist(13,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(2))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(2)))).^2;
Dist(14,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(3))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(2)))).^2;

Dist(15,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(2))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(3)))).^2;
Dist(16,2:end) = abs((Chan_Op(2:end)-QPSK_SYM(3))+(pred_coef*(Chan_Op(1:end-1)-QPSK_SYM(3)))).^2;
log_gamma = -Dist/(2*pred_var); % log gamma

 
log_gamma1 = log_gamma(:,1:num_bit); % branch metrices for component decoder 1
log_gamma2 = log_gamma(:,num_bit+1:end); % branch metrices for component decoder 2
 
% a priori LLR for component decoder 1 for 1st iteration
LLR = zeros(1,num_bit);

% iterative logMAP decoding
LLR = log_BCJR(LLR,log_gamma1,num_bit); % outputs extrinsic information
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %1

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %2

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %3

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %4

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %5

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %6

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR(LLR(intr_map),log_gamma2,num_bit); %7

LLR = log_BCJR(LLR(deintr_map),log_gamma1,num_bit);
LLR = log_BCJR_END(LLR(intr_map),log_gamma2,num_bit); % 8: outputs aposteriori probabilities

% hard decision 
LLR = LLR(deintr_map);
dec_data = LLR<0;

 % Calculating total bit errors
C_Ber = C_Ber + nnz(dec_data-a); 
end

BER = C_Ber/(sim_runs*num_bit)
toc()

