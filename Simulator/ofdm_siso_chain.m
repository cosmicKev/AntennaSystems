
function ofdm_siso_chain
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
N_OFDM_SYM=1e4;
m=4; %QPSK Mod.
samp_freq=15.36e6;
load pdp.mat
noise_variance = 1.*10.^(-EbN0./10)./(log2(m));


h_temp=zeros(NT*NR,Nc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:length(EbN0)    % EbN0_dB
    
    
    for k=1:N_OFDM_SYM   % number of OFDM symbols

    %################# Transmitter #################
        % data generation
        data=GenData(Nc_aval*2);
        
        %data modulation 
        tx_data= mod_data(data,m);
        
        %OFDM FRAMING
        tx_frame=zeros(1,1024);
        tx_frame(129:896)=tx_data;
        
             % fft modulation
             tx_time=ifft(tx_frame)*sqrt(Nc);
             
        %guard interval insertion
        tx_time_gi=[tx_time(Nc-Ng+1:Nc) tx_time];
        
    %################# Channel Model #################
         %multipath channel
        
        [ht11, hf11]=channel_gen(pdp,samp_freq, Nc);
         hf11_d=hf11(129:896);
         
        %awgn noise 
        noise=sqrt(noise_variance(p)/2)*(randn(1,Nc+Ng)+1i*randn(1,Nc+Ng));
        
        Rx_Signal=conv_s_h( tx_time_gi,ht11,pdp,Nc,samp_freq,tg)+noise;
        
        %received signal assuming a only a noise channel
        Rx_Signal_awgn=tx_time_gi+noise;
        
   %################# Receiver #################
   
        %GI removal + FFT
        Rx_Signal_gif=fft(Rx_Signal(Ng+1:end));
        Rx_data=Rx_Signal_gif(129:896);
        
        Rx_Signal_gif_awgn=fft(Rx_Signal_awgn(Ng+1:end));
        Rx_data_awgn=Rx_Signal_gif_awgn(129:896);
        
        % Equalization
        softdata=Rx_data.*conj(hf11_d)./abs(hf11_d).^2;
        
        %data demodulation
        harddata=demod_data(softdata,m, Nc_aval*2);
        harddata_awgn=demod_data(Rx_data_awgn,m, Nc_aval*2);
        
        %BER computation
        ber_ofdm(k)=compute_ber(harddata,data);
        ber_ofdm_awgn(k)=compute_ber(harddata_awgn,data);
    
    end  
    
    ber_each_eb(p)=sum(ber_ofdm)/N_OFDM_SYM
    ber_each_eb_awgn(p)=sum(ber_ofdm_awgn)/N_OFDM_SYM
      
end
  


%  Ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
semilogy(EbN0,ber_each_eb,'b*-','LineWidth',2);
hold on
semilogy(EbN0,ber_each_eb_awgn,'k-','LineWidth',2);
axis([0 32 10^-5 0.2])
grid on
xlabel('Eb/No, dB');
ylabel('BER');
legend('OFDM Fading Channel','OFDM AWGN')



%%%%%%%%%%%%%%% Auxiliar Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************** Data Generation *********************
    
function data=GenData(N_Data)
% Generates the N_Data length vector with the rand data that will included
data=round(rand(1,N_Data));


function data_symbol=mod_data(data, m)
% Data Modulation
%BPSK, m=2
%QPSK, m=4,
%16-QAM, m=16

switch m
    case 2
        data_symbol = -data*2+1; %mapping of 1->-1 and 0->1

    case 4
        % Coding of data bits in QPSK symbols - using Grey coding
        % (00->1+i; 01->1-i; 10->-1+i; 11->-1-i)
        % bit MS defines real polarity
        % bit LS defines imag polarity
        data_temp = reshape(data,2,length(data)/2);
        data_real = data_temp(1,:);
        data_imag = data_temp(2,:);
        data_symbol = sqrt(2)/2*((-1).^(data_real)+i*(-1).^(data_imag));
        
     case 16
        data_temp = reshape(data,4,length(data)/4);
        data_r1 = data_temp(1,:);
        data_i1 = data_temp(2,:);
        data_r2 = data_temp(3,:);
        data_i2 = data_temp(4,:);
        data_symbol = 2/sqrt(10).*(0.5*(-1).^(data_r2).*(-1).^(data_r1)+(-1).^( data_r1)+i.*(0.5*(-1).^(data_i2).*(-1).^(data_i1)+ (-1).^(data_i1)));
      otherwise
        helpdlg('Constellation size (m) not available');
       
end




%**************** Data De-modulation *********************

function decoded_data=demod_data(data_symbol, m, N_Data)

vect_IMAG = imag(data_symbol);
vect_REAL = real(data_symbol);
coder_type_value=0;
switch m
    case 2
        if (coder_type_value==0)
        	%hard decision
	      decoded_data= ceil(-vect_REAL./(1.0000000000001.*max(abs(vect_REAL))));
 
        else
        	%soft decision
	        decoded_data = vect_REAL; 
        end
        
    case 4
        % Decoding of data bits in QPSK symbols - using Grey coding
        % (1+i->00; 1-i->01; -1+i->10; -1-i->11)
        % real polarity defines bit MS
        % imag polarity defines bit LS
        if (coder_type_value==0)
        	%hard decision
	        vect_REAL_1 = ceil(-vect_REAL./(1.0000000000001.*max(abs(vect_REAL))));
        	vect_IMAG_1 = ceil(-vect_IMAG./(1.0000000000001.*max(abs(vect_IMAG))));
	        decoded_data = reshape([vect_REAL_1; vect_IMAG_1],1,N_Data);
        else
        	%soft decision
	        decoded_data = reshape([vect_REAL; vect_IMAG],1,N_Data);
        end
    case 16
        P_1= vect_REAL;
        P_2= vect_IMAG;
        P_3= abs(vect_REAL)-2/sqrt(10);
        P_4= abs(vect_IMAG)-2/sqrt(10);
        if (coder_type_value==0)
            %hard decision
            vect_IMAG_1 = ceil(-P_2./(1.0000000000001.*max(abs(P_2))));
            vect_IMAG_2 = ceil(-P_4./(1.0000000000001.*max(abs(P_4))));
            vect_REAL_1 = ceil(-P_1./(1.0000000000001.*max(abs(P_1))));
            vect_REAL_2 = ceil(-P_3./(1.0000000000001.*max(abs(P_3))));
            decoded_data =  reshape([vect_REAL_1; vect_IMAG_1; vect_REAL_2; vect_IMAG_2],1,N_Data);
        else    
            %soft decision
            decoded_data =  reshape([P_1; P_2; P_3; P_4],1,N_Data);
        end

    otherwise
        helpdlg('Constellation size (m) not available');
end
%*****************Channel Generation *****************
function [ht, Hf]=channel_gen(pdp,samp_freq, Nc)

delta_t=1/samp_freq;              %sample duratiion
Npaths = length(pdp(:,1));        % No. of paths considered for the channel
deltans=delta_t/1e-9;             % Sampling interval in ns

path_pot_lin=10.^(pdp(:,2)/10);
path_pot_lin=path_pot_lin./sum(path_pot_lin);

delays = pdp(:,1);
delays = round(delays./(deltans))+1;

multipath = zeros(1,Npaths);

for n=1:Npaths
pot_componente=0.5;
multipath(n)=sqrt(pot_componente)*randn(1,1)+j*sqrt(pot_componente)*randn(1,1);
multipath(n)=multipath(n).*sqrt(path_pot_lin(n));

end

RI=zeros(1,Nc);
RI(delays) = RI(delays) + multipath;

ht=RI;
Hf=fft(ht);

%****************CONV_S_H
function y=conv_s_h(s,h,pdp,Nc,samp_freq,tg)
%tg=5.21e-6;                       % guerad time -> 80 samples

delays=pdp(:,1);
delta_t=1/samp_freq;
Npaths = length(delays);                      % No. of paths considered for the channel
deltans=delta_t/1e-9;             % Sampling interval in ns


delays = round(delays/(deltans))+1;
Ng= round(tg/(delta_t))+1;
f_zeros=find(h==0);
h(f_zeros)=[];

aux = zeros(Npaths,Nc+Ng);

for n=1:Npaths
    

conv_sh=h(n)*s;

aux(n,delays(n):Nc+Ng-1+delays(n)-1)=conv_sh;

end


y=sum(aux);
y=y(1:Nc+Ng-1);

%**************** BIT error rate *********************

function [ber, Nerror]=compute_ber(x,y)

Nerror=sum(x~=y);
ber=Nerror/length(x);

