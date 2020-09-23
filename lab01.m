clear all
clc
prompt1 = 'Please enter modulation_scheme\n';
modulation = input(prompt1,'s');

prompt2 = 'Please enter snr: \n';
snr = input(prompt2);

prompt3 = 'Please enter number_of_packets: \n';
number_of_packets = input(prompt3);

prompt4 = 'Please enter bytes: \n';
bytes = input(prompt4);

run(modulation, snr, number_of_packets, bytes);

TX_POWER_SCALE = 1.0;
%snr = 30;
%modulation = '64QAM';
%packets = 100;
%bytes = 1440;

total_SER = 0;
total_syms = 0;
total_BER = 0;
total_bit = 0;
total_evm = 0;
storage_size = 0;
total_packetLoss = 0;

if strcmp(modulation, 'BPSK')
    storage_size = 8;
elseif strcmp(modulation, 'QPSK')
    storage_size = 4;
elseif strcmp(modulation, '16QAM') || strcmp(modulation, '64QAM')
    storage_size = 2;
else
    fprintf('Invaild input, must be either BPSK. QPSK, 16QAM or 64QAM');
    return;
end


% logics for dividing bytes to each packet
packetSize = zeros(number_of_packets,1);
for p = 1:number_of_packets
    if p == length(packetSize)
        if (mod(bytes, number_of_packets) ~= 1)
           packetSize(p) = floor(bytes/number_of_packets)+mod(bytes, number_of_packets);
        end
    else 
        packetSize(p) = floor(bytes/number_of_packets);
    end
end

tx_vec_real = zeros(1,1);
tx_vec_imag = zeros(1,1);
rx_vec_real = zeros(1,1);
rx_vec_imag = zeros(1,1);
index = 1;

for p = 1:length(packetSize)
    packet_SER = 0;
    packet_syms = 0;
    packet_BER = 0;
    packet_total_bit = 0;
    
    p_counter = 0;
    
    packet_Rx_syms_real = zeros(1,2*storage_size);
    packet_Rx_syms_imag = zeros(1,2*storage_size);
    packet_TX_syms_real = zeros(1,storage_size,1);
    packet_TX_syms_imag = zeros(1,storage_size,1);
    
    eachPacketSize = packetSize(p);
    
    % for sending a packet
    for a = 1:eachPacketSize
        randNum = randi(255);
        
        a_counter = 0;

        % convert from demical number to binary
        randBit = de2bi(randNum,'left-msb',8);
        if strcmp(modulation, 'BPSK')
            bitMap = reshape(randBit, 8, []);
        elseif strcmp(modulation, 'QPSK')
            bitMap = [randBit(1:2);randBit(3:4);randBit(5:6);randBit(7:8)];
        elseif strcmp(modulation, '16QAM')
            bitMap = [randBit(1:4);randBit(5:8)];
        elseif strcmp(modulation, '64QAM')
            bitMap = [[zeros(1,4),randBit(1:2)];randBit(3:8)];
        else
            fprintf('Invaild input, must be either BPSK. QPSK, 16QAM or 64QAM');
            return;
        end 
        tx_data = bi2de(bitMap,'left-msb');


        modvec_bpsk = (1/sqrt(2)).*[-1 1];
        modvec_16qam = (1/sqrt(10)).*[-3 -1 1 3];
        modvec_64qam = (1/sqrt(42)).*[-7 -5 -3 -1 1 3 5 7];

        mod_fcn_bpsk = @(x) complex(modvec_bpsk(1+x),0);
        mod_fcn_qpsk = @(x) complex(modvec_bpsk(1+bitshift(x,-1)), modvec_bpsk(1+mod(x,2)));
        mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x,-2)),modvec_16qam(1+mod(x,4)));
        mod_fcn_64qam = @(x) complex(modvec_64qam(1+bitshift(x,-3)), modvec_64qam(1+mod(x,8)));

        % Signal modulations for the data
        if strcmp(modulation, 'BPSK')
            tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
        elseif strcmp(modulation, 'QPSK')
            tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
        elseif strcmp(modulation, '16QAM')
            tx_syms = arrayfun(mod_fcn_16qam, tx_data);
        elseif strcmp(modulation, '64QAM')
            tx_syms = arrayfun(mod_fcn_64qam, tx_data);
        else
            fprintf('Invaild input, must be either BPSK. QPSK, 16QAM or 64QAM');
            return;
        end

        tx_length = length(tx_syms);
        tx_element_len = length(tx_syms(1,1:end));
        rx_syms = zeros(tx_length, 2*tx_element_len);

        for i = 1:tx_length
            tx_vec = ifft(tx_syms(i,1:end));
            %tx_vec_real= real(tx_vec);
            %tx_vec_imag= imag(tx_vec);
            
            tx_vec_real(index,1:end) = real(tx_vec);
            tx_vec_imag(index,1:end) = imag(tx_vec);            
            
            %disp(tx_data(i,1:end));
            %disp(tx_syms(i,1:end));

            % Define a half-band 2x interpolation filter response
            %interp_filt2 = zeros(1,43);
            %interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
            %interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
            %interp_filt2(22) = 16384;
            %interp_filt2 = interp_filt2./max(abs(interp_filt2));
            % interpolate
            %tx_vec_2x = zeros(1, 2*numel(tx_vec));
            %tx_vec_2x(1:2:end) = tx_vec;
            %tx_vec = filter(interp_filt2, 1, tx_vec_2x);
            %tx_vec = TX_POWER_SCALE.* tx_vec./ max(abs(tx_vec));

            % add Gaussian noise
            tx_channel = awgn(tx_vec,snr,'measured');

            %rx_syms_raw = filter(interp_filt2, 1, tx_channel);
            %rx_syms_raw = rx_syms_raw(1:2:end);

            % Receiver for BPSK
            rx_temp = fft(tx_channel);
            
            rx_vec_real(index,1:end) = real(rx_temp);
            rx_vec_imag(index,1:end) = imag(rx_temp);
            
            index = index+1;

            temp_len = length(rx_temp);
            for j = 1:temp_len
                rx_syms(i,j) = rx_temp(j);
            end
        end 

        a_counter = a_counter+1;
        demod_fcn_bpsk = @(x) double(real(x)>0);
        demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
        demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));
        demod_fcn_64qam = @(x) (32*(real(x)>0)) + (16*(abs(real(x))<0.6172)) + (8*((abs(real(x))<(0.9258))&&((abs(real(x))>(0.3086))))) + (4*(imag(x)>0)) + (2*(abs(imag(x))<0.6172)) + (1*((abs(imag(x))<(0.9258))&&((abs(imag(x))>(0.3086)))));

        % Signal demodulation depends on the modulation scheme
        if strcmp(modulation, 'BPSK')
            rx_data = arrayfun(demod_fcn_bpsk, rx_syms);
        elseif strcmp(modulation, 'QPSK')
            rx_data = arrayfun(demod_fcn_qpsk, rx_syms);
        elseif strcmp(modulation, '16QAM')
            rx_data = arrayfun(demod_fcn_16qam, rx_syms);
        elseif strcmp(modulation, '64QAM')
            rx_data = arrayfun(demod_fcn_64qam, rx_syms);
        else
            fprintf('Invaild input, must be either BPSK. QPSK, 16QAM or 64QAM');
            return;
        end

        %disp(rx_syms);

        % Plots the TX and RX symbols
        rxSyms_real = real(reshape(rx_syms,1,[]));
        rxSyms_imag = imag(reshape(rx_syms,1,[]));
        txSyms_real = real(reshape(tx_syms,1,[]));
        txSyms_imag = imag(reshape(tx_syms,1,[]));

        packet_Rx_syms_real(a,1:end) = rxSyms_real;
        packet_Rx_syms_imag(a,1:end) = rxSyms_imag;
        packet_TX_syms_real(a,1:end) = txSyms_real;
        packet_TX_syms_imag(a,1:end) = txSyms_imag;
        % plot(rxSyms_real,rxSyms_imag,'r.',txSyms_real,txSyms_imag,'bo');
        %grid
        %legend('RX','TX')
        %title('TX and RX Constellations')

        % Calculations for SER and total SER
        ser = sum(tx_data ~= rx_data(:,1));
        packet_SER = packet_SER+ser;
        sizeRx = size(rx_syms);
        packet_syms = packet_syms+sizeRx(1);

        % Calculations for BER and total BER
        decBin = bitxor(de2bi(tx_data,'left-msb',8), de2bi(rx_data(:,1),'left-msb',8));
        bitXor = find(reshape(decBin,1,[]));
        ber = size(bitXor);
        packet_BER = packet_BER+ber(2);
        packet_total_bit = packet_total_bit+8;

        %fprintf('Sym Errors:  %d \n', ser);
        %fprintf('Bit Errors:  %d \n', ber(2));
    end 
    
    p_counter = p_counter+1;
    
    total_SER = total_SER+packet_SER;
    total_syms = total_syms+packet_syms;
    total_BER = total_BER+packet_BER;
    total_bit = total_bit+packet_total_bit; 
    packet_evm = sqrt((sum((packet_Rx_syms_real(:,(1:storage_size))-packet_TX_syms_real).^2+(packet_Rx_syms_imag(:,(1:storage_size))-packet_TX_syms_imag)).^2)/packet_syms);
    total_evm = total_evm + packet_evm;
    
    if (packet_SER > 0)
        total_packetLoss = total_packetLoss+1;
    end
    
    plot(packet_Rx_syms_real(:,(1:storage_size)),packet_Rx_syms_imag(:,(1:storage_size)),'r.',packet_TX_syms_real,packet_TX_syms_imag,'bo');
    title('TX and RX Constellations')
    if (p == 1)
        hold on  
    elseif (p == length(packetSize))
        hold off 
        grid
        legend('Rx','Tx');
        figure
    end
end

subplot(2,1,1);
plot((1:length(tx_vec_real)), tx_vec_real);
title('In-phase TX waveform');
subplot(2,1,2); 
plot((1:length(tx_vec_imag)), tx_vec_imag);
title('Quadrature TX waveform');
figure 

subplot(2,1,1);
plot((1:length(rx_vec_real)), rx_vec_real);
title('In-phase RX waveform');
subplot(2,1,2); 
plot((1:length(rx_vec_imag)), rx_vec_imag);
title('Quadrature RX waveform');


fprintf('Total Symbol Errors:  %f \n', (total_SER(1)/total_syms));
fprintf('Total Bit Errors:  %f \n', (total_BER/total_bit));
fprintf('RX EVM: %f\n',sum(total_evm/total_syms)/storage_size);
fprintf('Packet loss: %f \n',total_packetLoss/number_of_packets);


