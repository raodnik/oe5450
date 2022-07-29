clc;
clear;

V = readmatrix('p.txt');

imagesc(V)
set(gca, 'YTickLabel', 0.9:-0.1:0)
set(gca, 'XTickLabel', 0.1:0.1:1)