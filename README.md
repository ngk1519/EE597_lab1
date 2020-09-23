# EE597_lab1

Author: Kevin Hock Tuen Ng, Xinhong Liu

## Project description 

This project was one of the lab assignments from the course EE597 Wireless Networks at the University of Southern California, Viterbi School of Engineering. 

The goal was to simulate and plot digital modulations for wireless communications with various modulation schemes, and calculate the bit and symbol error rates. The user can select one of the following schemes: Binary Phase-shift keying (BPSK), Quadrature Phase-shift keying (QPSK), 16 Quadrature amplitude modulation (16QAM) and 64 Quadrature amplitude modulation (64QAM). 

## Source files

1. lab01.m
	This file produces the constellation diagram, transmitter and receiver waveforms based on the user's selections of modulation scheme, noise SNR, number of packets and total transmit byte size. 

2. run.m
	This file allows users to input desired settings for the simulation from the lab01.m file.

3. lab1_v3.m
	As this file simulates 4 different modulation schemes (BPSK, QPSK, 16QAM, 64QAM), file computes the bit error rate (BER), symbol error rate (SER), error vector magnitude (EVM) and packet loss percentages for various levels of noise power in the unit of dBm. 

## References/Credits

1. WARP example: https://warpproject.org/trac/browser/ResearchApps/PHY/WARPLAB/WARPLab7/M_Code_Examples/wl_example_siso_ofdm_txrx.m
