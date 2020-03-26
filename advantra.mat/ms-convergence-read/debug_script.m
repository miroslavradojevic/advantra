clear all; close all; clc;
% set these!
a = csvread('~/out.csv');
im = imread('~/test2.tif');
iter_nr = 300;
point = [30 8];



h = size(im, 1);
w = size(im, 2);
img_size = h*w;

figure;
imshow(im); hold on;
plot(point(2), point(1), 'go', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g','MarkerSize',5);
plot(a(1:iter_nr,4),                a(1:iter_nr,3),             'r-');
plot(a(iter_nr+1:2*iter_nr,4),      a(iter_nr+1:2*iter_nr,3),   'gx');
plot(a(2*iter_nr+1:3*iter_nr,4),    a(2*iter_nr+1:3*iter_nr,3), 'b.');
%legend('start','rnd. before sampling values', 'blurring', 'rnd. after sampling values');

start = 3*iter_nr;
figure;
imshow(im); hold on;
plot(a(start+1:    start+img_size,4),    a(start+1:     start+img_size,3), 'r.');

start = start+img_size;
figure;
imshow(im); hold on;
plot(a(start+1:    start+img_size,4),    a(start+1:     start+img_size,3), 'g.');

start = start+img_size;
figure;
imshow(im); hold on;
plot(a(start+1:    start+img_size,4),    a(start+1:     start+img_size,3), 'b.');

figure;
plot(1:iter_nr, a(1:iter_nr, 2), 'r'); grid on; hold on;
plot(1:iter_nr, a(iter_nr+1:2*iter_nr, 2), 'g');
plot(1:iter_nr, a(2*iter_nr+1:3*iter_nr, 2), 'b');

disp('done');