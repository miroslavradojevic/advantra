clear all; close all; clc;

a1 = csvread('~/out8.csv');
im = imread('~/test2.tif');

figure;
imshow(im); hold on;
plot(a1(:,2), a1(:,1), 'rx');
export_fig('~/out8.png', '-png', '-a1', '-native');
% figure;
% imshow(im); hold on;
% plot(a2(:,2), a2(:,1), 'rx');
% export_fig('~/2.png', '-png', '-a1', '-native');
% figure;
% imshow(im); hold on;
% plot(a3(:,2), a3(:,1), 'rx');
% export_fig('~/3.png', '-png', '-a1', '-native');
% figure;
% imshow(im); hold on;
% plot(a4(:,2), a4(:,1), 'rx');
% export_fig('~/4.png', '-png', '-a1', '-native');