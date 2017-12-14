function [sn] = sf_decoding(rx_deframing_data,h11,h21)
    yn = rx_deframing_data;
    sn = zeros(1,768);
    sn(1:2:end) = 1/sqrt(2) * conj(h11(2:2:end)).*yn(1:2:end) + 1/sqrt(2) * h21(1:2:end).*conj(yn(2:2:end));
    sn(2:2:end) = 1/sqrt(2) * conj(h21(2:2:end)).*yn(1:2:end) - 1/sqrt(2) * h11(1:2:end).*conj(yn(2:2:end));
end