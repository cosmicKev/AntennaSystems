function [ tx_time_gi ] = ofdmFraming( sfcTx, Nc, Ng )
tx_frame=zeros(1,1024);
tx_frame(129:896)=sfcTx;

% fft modulation
tx_time=ifft(tx_frame)*sqrt(Nc);
%guard interval insertion
tx_time_gi=[tx_time(Nc-Ng+1:Nc) tx_time];

end

