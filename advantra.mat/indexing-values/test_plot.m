close all; clear all; clc;

%		double step 	= 0.01;// how many test samples to take into account
%		double start 	= Math.PI;		
%		double end 	= 2*Math.PI;
%		int N  = 10; // levels


%%%%%%%% works with vals_incl, vals_excl (taken from Java) %%%%%%%%%%

vals; % to load them

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure(1); 
title('last el. INcluded');
plot(vals_incl(:,2), vals_incl(:,1), 'ro'); hold on; 
plot(vals_incl(:,2), vals_incl(:,3), 'bx'); grid on;
xlabel('index'); 
ylabel('index covered values');
legend('original', 'corr. discrete(N=10)');

print -deps -color included-index-coverage.eps
close;

figure; 
title('last interval el. INcluded');
plot(1:size(vals_incl,1), vals_incl(:,1), 'ro'); hold on; 
plot(1:size(vals_incl,1), vals_incl(:,3), 'bx'); grid on;
xlabel(''); 
ylabel('values');
legend('original', 'corr. discrete(N=10)');

print -deps -color included-just-values.eps
close;

figure; 
title('last el. EXcluded');
plot(vals_excl(:,2), vals_excl(:,1), 'ro'); hold on; 
plot(vals_excl(:,2), vals_excl(:,3), 'bx'); grid on;
xlabel('index'); 
ylabel('index covered values');
legend('original', 'corr. discrete(N=10)');

print -deps -color excluded-index-coverage.eps
close;

figure; 
title('last interval el. INcluded');
plot(1:size(vals_excl,1), vals_excl(:,1), 'ro'); hold on; 
plot(1:size(vals_excl,1), vals_excl(:,3), 'bx'); grid on;
xlabel(''); 
ylabel('values');
legend('original', 'corr. discrete(N=10)');

print -deps -color excluded-just-values.eps
close all; 

disp('done. plots saved.');
