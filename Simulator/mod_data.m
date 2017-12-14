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