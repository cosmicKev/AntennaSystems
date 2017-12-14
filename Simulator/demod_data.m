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
