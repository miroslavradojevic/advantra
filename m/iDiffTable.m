clear all; close all; clc;

disp("SNR\tIDIFF");

bg = 30;

for snr = 1 : 5,

	idiff = (snr + sqrt(snr^2+4*snr*bg)) / 2;
	disp(strcat(num2str(snr), " \t " , num2str(idiff)));

end
