% ---- Randon Seed ----
s = rng(311);       % Seed the RNG for repeatability

% ---- Code parameters ----
K = 54;             % Message length in bits, including CRC, K > 24 for DL
E = 124;            % Rate matched output length, E <= 8192

EbNo = 0.8;         % EbNo in dB
L = 8;              % List length, a power of two, [2 4 8]
numFrames = 10;     % Number of frames to simulate

% ---- Code Construction ----
% Downlink channel parameters (K > 24, including CRC bits)
crcLen = 24;      % Number of CRC bits for DL, Section 5.1, [7]
nPC = 0;          % Number of parity check bits, Section 5.3.1.2, [7]
nMax = 9;         % Maximum value of n, for 2^n, Section 7.3.3, [7]
iIL = true;       % Interleave input, Section 5.3.1.1, [7]
iBIL = false;     % Interleave coded bits, Section 5.4.1.3, [7]

% Code construction
F = h5gPolarConstruct(K,E,nMax);  % 0 for information, 1 for frozen
N = length(F);                    % Mother code block length

R = K/E;                          % Effective code rate
bps = 2;                          % bits per symbol, 1 for BPSK, 2 for QPSK
EsNo = EbNo + 10*log10(bps);
snrdB = EsNo + 10*log10(R);       % in dB
noiseVar = 1./(10.^(snrdB/10));

% ---- Polar Encoding ----
% Polar Encoder
polarEnc = h5gPolarEncoder(N,K,F,'InterleaveInput',iIL);

% Modulator, Channel, Demodulator
qpskMod = comm.QPSKModulator('BitInput', true);
chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);
qpskDemod = comm.QPSKDemodulator('BitOutput',true,'DecisionMethod', ...
    'Approximate log-likelihood ratio','Variance',noiseVar);

% Polar Decoder
polarDec = h5gPolarDecoder(N,K,F,L,crcLen,'DeinterleaveOutput',iIL);

% Bit-Error rate meter
ber = comm.ErrorRate;

numferr = 0;
for i = 1:numFrames

    % Generate a random message
    msg = randi([0 1],K-crcLen,1);

    % CRC attachment
    msgcrc = h5gCRCEncode(msg,'24C');

    % Polar encode
    encOut = polarEnc(msgcrc);

    % Rate match
    modIn = h5gRateMatchPolar(encOut,K,E,iBIL);

    % Modulate
    modOut = qpskMod(modIn);

    % Add white Gaussian noise
    rSig = chan(modOut);

    % Soft demodulate
    rxLLR = qpskDemod(rSig);

    % Rate recover
    decIn = h5gRateRecoverPolar(rxLLR,K,N,iBIL);

    % Polar decode
    decBits = polarDec(decIn);

    % Compare msg and decoded bits
    errStats = ber(decBits(1:K-crcLen), msg);
    numferr = numferr + any(decBits(1:K-crcLen)~=msg);

end

disp(['Block Error Rate: ' num2str(numferr/numFrames) ...
      ', Bit Error Rate: ' num2str(errStats(1)) ...
      ', at SNR = ' num2str(snrdB) ' dB'])

rng(s);     % Restore RNG
