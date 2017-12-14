function [ Rx_data ] = ofdmDeframing( Rx_Signal, Ng )
    %GI removal + FFT
        Rx_Signal_gif=fft(Rx_Signal(Ng+1:end));
        Rx_data=Rx_Signal_gif(129:896);
end

