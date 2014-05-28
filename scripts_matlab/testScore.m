close all; clear all; clc;

disp('start...');

nd = 3;
sc = 1.0;
r1 = nd*sc;
r2 = nd*(sc+0.5);
r3 = nd*(sc+0.5+0.5);
dCirc = 2.0;

disp('generating points...');
tripletsX = [];
tripletsY = [];
cnt = 0;
for o1 = 0 : dCirc : 2*pi*r1-dCirc,
	%o1 = o1 + dCirc;
	a1 = o1/r1;
	%disp(a1);
	for o2 = 0 : dCirc : 2*pi*r2-dCirc,
		%o2 = o2 + dCirc;
		a2 = o2/r2;
			for o3 = 0 : dCirc : 2*pi*r3-dCirc,
				%o3 = o3 + dCirc;
				a3 = o3/r3;
				cnt = cnt + 1;	
				tripletsX = [tripletsX; r1*cos(a1), r2*cos(a2), r3*cos(a3)];
				tripletsY = [tripletsY; r1*sin(a1), r2*sin(a2), r3*sin(a3)];
			end
	end
end
disp('total:');
disp(cnt);

%figure;
%plot(0, 0, 'bo'); grid on; axis equal; hold on;
%for i = 1 : 1 : size(tripletsX, 1),
%	plot(tripletsX(i,:), tripletsY(i,:), 'ro');
%end

disp(' calculate scores.. ');

A = zeros(size(tripletsX, 1), 1);
B = zeros(size(tripletsX, 1), 1);

for i = 1 : 1 : size(tripletsX, 1),
	x21 = tripletsX(i, 2) - tripletsX(i, 1);
	x32 = tripletsX(i, 3) - tripletsX(i, 2);
	y21 = tripletsY(i, 2) - tripletsY(i, 1);		
	y32 = tripletsY(i, 3) - tripletsY(i, 2);	
	A(i) = (x21*x32+y21*y32)/(sqrt(x21^2+y21^2)*sqrt(x32^2+y32^2)); % cos
	A(i) = 0.5*(A(i)+1);
	x10 = tripletsX(i, 1) - 0;
	x31 = tripletsX(i, 3) - tripletsX(i, 1);
	y10 = tripletsY(i, 1) - 0;
	y31 = tripletsY(i, 3) - tripletsY(i, 1);
	B(i) = (x10*x31+y10*y31)/(sqrt(x10^2+y10^2)*sqrt(x31^2+y31^2)); % cos
	B(i) = 0.5*(B(i)+1);	
end

C = A;
%C = A.*B;

[Csorted, Isorted] = sort(C, 'descend'); %descend


figure;
plot(0, 0, 'bo'); grid on; axis equal; hold on;
for i = 10 : 1 : 20,
	plot(tripletsX(Isorted(i),:), tripletsY(Isorted(i),:), 'r-');
end


