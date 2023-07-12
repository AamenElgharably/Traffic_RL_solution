% LTE Physical Layer simulation with AWGN channel, Turbo coding according to 20MHz 4G specifications
%% Initialization
   clear all; close all; clc;
%%  Define system parameters
%numSubcarriers =1200;           % Total Number of subcarriers 
NFFT=2048;                      %FFT size
CP = 512;                       % Long Cyclic prefix 
CP_start =NFFT -CP ;            % Cyclic prefix Start
CP_end =  NFFT;                 %Cyclic prefix end
df = 15e3;                      % Subcarrier spacing
fs = NFFT*df;                   % Sampling frequency
Ts = 1/fs;                      % Sampling period
Total_Length = NFFT+CP;         % OFDM symbol size
numSymbols = 7;                 % Number of OFDM symbols
numSlots = 2;                   % number of slots
numSymbolSlot = 7;              % number of symbols per slot
numsubCperRB = 12;              % number of subcarriers per resource block
QAMorderList = [4,16,64];       % the value of M in M-QAM, modulation used for data
M=64;                            
codeRate = 1/3;                 %The code rate for Turbo coding
K = log2(M);                    % Number of bits per symbol
%% Tx
load('UE2.mat');
UE=UE2;
%% Encoding, and Decoding

UE_nums=length(UE);
FRM=[128 256 512 1024 2048 4096];
for i=1:UE_nums
    
    maxBits = UE(i).Allocated_RBs*numsubCperRB*numSymbolSlot*numSlots;
    for k=1:length(FRM)
        if (FRM(k)*3)/UE(i).Modulation < maxBits
            UE(i).FRM= FRM(k);
        end
    end    
end
for i=1:UE_nums
FRM = UE(i).FRM;                     %Frame size, including the CRC bits 128,256,512,1024,2048,6144 
FRM = FRM-24;                   %The frame size, excluding the CRC bits.
K1 = FRM+24;                    %The length of the interleaver.
[f1, f2] = getf1f2(K1);        %f1, f2: Parameters used to generate the interleaver indices 
inx=0:(K1-1);   
Indices = mod(f1*inx + f2*inx.^2, K1) + 1; %indices: The interleaver indices used for Turbo coding.
Trellis = poly2trellis(4, [13 15], 13); %The trellis structure used for Turbo coding.
%{
n_FRM: This parameter specifies the number of frames to transmit. Each frame contains a block of random data and is processed independently.

FRM: This parameter specifies the size of each frame, excluding the CRC bits. The value of FRM is set to 2048, which is the standard frame size used in the LTE system.

FRM = FRM-24: This line adjusts the value of FRM to include the 24-bit CRC (Cyclic Redundancy Check) that is added to each frame to enable error detection.

K1: This parameter specifies the length of the interleaver used in the Turbo coding process. The value of K1 is set to FRM+24, which is the size of the frame including the CRC bits.

codeRate: This parameter specifies the code rate used in the Turbo coding process. The value of codeRate is set to 1/3, which means that each input bit is encoded into three output bits.

f1=263; f2=480; inx=0:(K1-1);: These lines set the parameters used to generate the interleaver indices for the Turbo coding process. The variable inx is a vector of integers from 0 to K1-1.

Indices = rem((f1*inx)+(f2*(inx.^2)),K1)+1;: This line generates the interleaver indices used in the Turbo coding process. The indices are calculated using a simple formula that involves the values of f1, f2, and inx.

Trellis = poly2trellis(4, [13 15], 13);: This line defines the trellis structure used in the Turbo coding process. The poly2trellis function takes two arguments: the constraint length of the trellis (4 in this case), and the generator polynomials for the constituent convolutional codes ([13 15] in this case). The last argument specifies the starting state of the encoder.

%}
% Turbo encoder object
UE(i).TurboEnc = comm.TurboEncoder('TrellisStructure',Trellis,'InterleaverIndices',Indices);
% Turbo decoder object
UE(i).TurboDec = comm.TurboDecoder('TrellisStructure',Trellis,'InterleaverIndices',Indices,'NumIterations',3);
%{
Together, the Turbo encoder and decoder objects implement the Turbo coding scheme used in the LTE physical layer simulation.
The encoder object takes in a vector of input bits, adds redundant parity bits using the Turbo coding scheme,
and produces a vector of encoded bits. The decoder object takes in a vector of received bits,
applies the Turbo decoding process to correct any errors in the received signal, and produces a vector of decoded bits.

The TrellisStructure and InterleaverIndices parameters used in the encoder and decoder objects must be the same to ensure that the encoding 
and decoding processes are consistent.
The NumIterations parameter used in the decoder object specifies the number of iterations used in the decoding process and can affect the decoding performance.
In this code, 6 iterations are used for the decoding process.
%}    
end
CRC_Tx = comm.CRCGenerator('Polynomial',[1 1 zeros(1, 16) 1 1 0 0 0 1 1]);
% CRC detector object
CRC_Rx = comm.CRCDetector('Polynomial', [1 1 zeros(1, 16) 1 1 0 0 0 1 1]);
%{
CRC_Tx = comm.CRCGenerator('Polynomial',[1 1 zeros(1, 16) 1 1 0 0 0 1 1]);: This line creates a CRC generator object that is used to compute the CRC bits for each frame of data. The Polynomial parameter specifies the generator polynomial used to generate the CRC bits.

CRC_Rx = comm.CRCDetector('Polynomial', [1 1 zeros(1, 16) 1 1 0 0 0 1 1]);: This line creates a CRC detector object that is used to check the CRC bits received in each frame and detect any errors. The Polynomial parameter must match the generator polynomial used in the CRC generator object.
%}
%% Defining Resource Grid (RG) parameters
numRB = 100;                             % number of resource blocks allocated 
numsubC = numRB*numsubCperRB;           % number of sbcarriers available for transmission 
numSymbols = numSlots*numSymbolSlot;    % Total number of OFDM symbols in the slot
% Generating the resource grid 
RG = zeros(numsubC, numSymbols);

%% Generating Random data to test the previous part
  numBits = numsubC*K*numSymbolSlot*numSlots; % Computing the number of bits necessary
  dataBits = randi([0,1], numBits, 1);  % Generating the bits
  modulatedSymbols = qammod(dataBits, M, 'gray' , "InputType","bit", "UnitAveragePower", true);   % Generating QAM symbols 
  RG = reshape(modulatedSymbols, size(RG));   % Generating QAM symbols 
  plotResourceGrid(abs(RG), "RG", "OFDM symbols", "Subcarriers");  % Plotting the generated resource grid

%% Signal Source 
for i=1:UE_nums
  UE(i).Data= randi([0 1], UE(i).FRM-24,1); % Random binary data stream
end
%% CRC 
%t_data_crc = step(CRC_Tx,t_data);
for i=1:UE_nums
    UEPL(i).t_data_crc=step(CRC_Tx,UE(i).Data);
end    
%% Turbo Coding_Total_Length
%cod_data=step(TurboEnc,t_data_crc);
for i=1:UE_nums
    UEPL(i).cod_data=step(UE(i).TurboEnc, UEPL(i).t_data_crc);
end    
%% Rate Matching
D = K1 +4;
for i=1:UE_nums
  if UE(i).R==(UE(i).R);
    UEPL(i).ek_Length=length(UEPL(i).cod_data); 
  else
    UEPL(i).ek_Length=UE(i).Modulation*ceil((D/UE(i).Modulation)/UE(i).R);
  end 
UEPL(i).cod_data_matched=UEPL(i).cod_data(1:UEPL(i).ek_Length);
end
%% Modulation
%mod_data = qammod(cod_data_matched, M, 'gray' , "InputType","bit", "UnitAveragePower", true);   % Generating QAM symbols 
for i=1:UE_nums
   UEPL(i).mod_data = qammod(UEPL(i).cod_data_matched, 2^UE(i).Modulation, 'gray' , "InputType","bit", "UnitAveragePower", true);   % Generating QAM symbols 
end
%% OFDMA  
for i = 1:UE_nums
       no_OFDM_Symbol=14;
       UEPL(i).Allocated_RBs=UE(i).Allocated_RBs;
       UEPL(i).numsubC=UEPL(i).Allocated_RBs*numsubCperRB;

   
    UEPL(i).numZeros = no_OFDM_Symbol*UEPL(i).numsubC - length(UEPL(i).mod_data);
      if UEPL(i).numZeros > 0
        UEPL(i).mod_data = [UEPL(i).mod_data; zeros(UEPL(i).numZeros, 1)]; % Pad mod_data with zeros
      end
    UEPL(i).mod_data_reshaped = reshape(UEPL(i).mod_data, UEPL(i).numsubC,no_OFDM_Symbol); % serial to parallel conversion 
  
end
%%
for i = 1:UE_nums
    %A(i)=UEPL(i).no_OFDM_Symbol;
    B(i)=UE(i).Priority;
end    
numSlots=2;
RG = zeros(numRB*numsubCperRB, numSymbolSlot*numSlots); %% Create a matrix of zeros with dimensions numsubC x numSymbols to represent the resource grid
%Fill the resource grid with the reshaped modulated data
space_Sc=0;
space_Sy=0;
subFRM=1;
%while (U_Allocated<length(A))
for k = 1:9
    if  space_Sc < numRB*numsubCperRB
        index = find(B == k);
          for l=1:length(index)
               RG(space_Sc+1:size(UEPL(index(l)).mod_data_reshaped,1)+space_Sc, space_Sy+1:(size(UEPL(index(l)).mod_data_reshaped,2))*subFRM) = UEPL(index(l)).mod_data_reshaped;
               space_Sc=space_Sc+size(UEPL(index(l)).mod_data_reshaped,1);
               if  space_Sc+228 >= 0.99*numRB*numsubCperRB                 
                 space_Sy=14;  
                 space_Sc=0;
                 subFRM=subFRM+1;
               end    
          end 
    end
end
plotResourceGrid(abs(RG), "RG", "OFDM symbols", "Subcarriers");% Plot the resource grid

%%
 % IFFT
 for i = 1:UE_nums 
   FFTgrid = zeros(NFFT, numSymbolSlot*numSlots);
   RG = zeros(numRB*numsubCperRB, numSymbolSlot*numSlots);
   RG(1:size(UEPL(i).mod_data_reshaped,1), 1:size(UEPL(i).mod_data_reshaped,2)) = UEPL(i).mod_data_reshaped;
   orgIndexSet = 0:numRB*numsubCperRB-1;
   newIndexSet = mod(orgIndexSet - numRB*numsubCperRB/2, NFFT) + 1;
   FFTgrid(newIndexSet,:) = RG;
   
  UEPL(i).IFFT_data_parallel = ifft(FFTgrid , NFFT, 1); 
 end  
  % parallel to serial coversion & addition cyclic prefix
    
for i = 1:UE_nums
    UEPL(i).IFFT_data_Serial = [];
    IFFT_data_Symbol=[];
    for k = 1:no_OFDM_Symbol
        IFFT_data_Symbol = UEPL(i).IFFT_data_parallel(:,k);
        UEPL(i).IFFT_data_Serial = [UEPL(i).IFFT_data_Serial; IFFT_data_Symbol(end-CP + 1:end,:); IFFT_data_Symbol];
    end
end
%% Channel
% Adding AWGN noise
SNR = 20; % Signal to noise ratio (dB)
for i = 1:UE_nums
%{
signalPower = mean(mean(abs(UEPL(i).IFFT_data_Serial).^2));
noisePower = signalPower / (10^(SNR/10));
noise = sqrt(noisePower/2) * (randn(size(UEPL(i).IFFT_data_Serial)) + 1i*randn(size(UEPL(i).IFFT_data_Serial)));
UEPL(i).channelOutput = UEPL(i).IFFT_data_Serial + noise;
%}
UEPL(i).channelOutput = awgn(UEPL(i).IFFT_data_Serial,SNR,'measured');
end
%% Rx.
%% serial to parallel coversion & Removing cyclic prefix
for i = 1:UE_nums
 CP_cumulative = 0;
 UEPL(i).FFT_data_parallel = zeros(NFFT,no_OFDM_Symbol);
 for k = 1:no_OFDM_Symbol
     CP_cumulative = CP_cumulative + CP;
     RXtimeDataSymbol = UEPL(i).channelOutput(CP_cumulative + (k-1)*NFFT + 1: CP_cumulative +k*NFFT,:);
     UEPL(i).FFT_data_parallel(:,k) = RXtimeDataSymbol;
 end 
end
%% FFT
for i = 1:UE_nums
 UEPL(i).mod_data_rx_reshaped = fft(UEPL(i).FFT_data_parallel, NFFT, 1);
end 
%% Parallel to Serial
    orgIndexSet = 0:numRB*numsubCperRB-1;
    newIndexSet = mod(orgIndexSet - numRB*numsubCperRB/2, NFFT) + 1;
for i = 1:UE_nums
  RXfreqData = UEPL(i).mod_data_rx_reshaped(newIndexSet, :);
  UEPL(i).mod_data_rx = RXfreqData(1:size(UEPL(i).mod_data_reshaped,1), :) ;
  UEPL(i).mod_data_rx_Serial = reshape(UEPL(i).mod_data_rx, [],1);
  UEPL(i).mod_data_rx_Serial= UEPL(i).mod_data_rx_Serial(1:end-UEPL(i).numZeros ,1);
end 
%% Demodulation
for i = 1:UE_nums
  UEPL(i).demod_Data = qamdemod(UEPL(i).mod_data_rx_Serial,2^UE(i).Modulation, 'gray', "OutputType","bit", "UnitAveragePower", true); 
end   
%% Rate de-matching
for i = 1:UE_nums
UEPL(i).demod_Data_Dematched = zeros(length(UEPL(i).cod_data),1);
UEPL(i).demod_Data_Dematched(1:UEPL(i).ek_Length) = UEPL(i).demod_Data;
  if UE(i).R == (UE(i).R)
    UEPL(i).demod_Data_Dematched = UEPL(i).demod_Data_Dematched(1:end);
  else
    UEPL(i).demod_Data_Dematched = UEPL(i).demod_Data_Dematched(1:UE(i).Modulation*(D/UE(i).Modulation));
  end
end  
%% Turbo - Decoding
for i = 1:UE_nums
  UEPL(i).R_data=step(UE(i).TurboDec,UEPL(i).demod_Data_Dematched);
end  
%% CRC
for i = 1:UE_nums
  UEPL(i).R_data_crc = step(CRC_Rx,UEPL(i).R_data);
end  
%% BER Computation
for i = 1:UE_nums
%[number_of_errors(i),bit_error_rate_temp(i)] = biterr(UE(i).Data,UEPL(i).R_data_crc) ;
 UEPL(i).BER = sum(UE(i).Data ~= UEPL(i).R_data_crc)/length(UE(i).Data);
end 
   