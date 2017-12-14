function trab_mimo_ofdm_67555
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference ofdm chain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
% allocating memory & Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NT=2;  % number of TX antennas 
NR=1;  % number of RX antennas
Nc=1024; % number of subcarriers
Nc_aval=768;
Ng=80;   %guard interval
tg=5.21e-6;
EbN0=0:4:32;
%EbN0=100
N_OFDM_SYM=1e3;
m=4; %QPSK Mod.
samp_freq=15.36e6;
load pdp.mat
noise_variance = 1.*10.^(-EbN0./10)./(log2(m));

h_temp=zeros(NT*NR,Nc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:length(EbN0)    % EbN0_dB
    for k=1:N_OFDM_SYM   % number of OFDM symbols
    %################# Transmitter #################
        % Data Generation
            data=GenData(Nc_aval*2);
        % Data Modulation
            tx_data= mod_data(data,m);
        % Data Space-Frequency Coding
            [sfcTx1, sfcTx2] = sf_coding(tx_data);
        % OFDM FRAMING 1
            tx1_time_gi = ofdmFraming(sfcTx1, Nc, Ng);
        % OFDM FRAMING 2
            tx2_time_gi = ofdmFraming(sfcTx2, Nc, Ng);
        
    %################# Channel Model #################
        % Channel
            [ht11, hf11]=channel_gen(pdp,samp_freq, Nc);
             hf11_d=hf11(129:896);
            [ht21, hf21]=channel_gen(pdp,samp_freq, Nc);
             hf21_d=hf21(129:896);
            % AWGN Noise
                noise=sqrt(noise_variance(p)/2)*(randn(1,Nc+Ng)+1i*randn(1,Nc+Ng));
            Rx_Signal = conv_s_h(tx1_time_gi,ht11,pdp,Nc,samp_freq,tg) + conv_s_h(tx2_time_gi,ht21,pdp,Nc,samp_freq,tg) + noise;
        
   %################# Receiver #################
        % Single Antenna Receiver
        % OFDM_FRAMING 
            Rx_data = ofdmDeframing(Rx_Signal, Ng);
        % Data Space-Frequency Decoding
            softdata = sf_decoding(Rx_data,hf11_d,hf21_d);
        % data demodulation
             harddata=demod_data(softdata,m, Nc_aval*2);
        %BER computation
            ber_ofdm(k)=compute_ber(harddata,data);
         % ber_ofdm_awgn(k)=compute_ber(harddata_awgn,data);
    end
     ber_each_eb(p)=sum(ber_ofdm)/N_OFDM_SYM;
    % ber_each_eb_awgn(p)=sum(ber_ofdm_awgn)/N_OFDM_SYM
      
end
  
%  Ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
semilogy(EbN0,ber_each_eb,'b*-','LineWidth',2);
%hold on
%semilogy(EbN0,ber_each_eb_awgn,'k-','LineWidth',2);
axis([0 32 10^-5 0.2])
grid on
xlabel('Eb/No, dB');
ylabel('BER');
legend('OFDM Fading Channel','OFDM AWGN')

%%%%%%%%%%%%%%% Auxiliar Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    
   



%****************CONV_S_H


%**************** BIT error rate *********************
function [ber, Nerror]=compute_ber(x,y)

Nerror=sum(x~=y);
ber=Nerror/length(x);

