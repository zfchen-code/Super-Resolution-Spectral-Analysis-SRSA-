clc; clear; close all;

rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir,'examples'));
addpath(fullfile(rootDir,'src_SRSA'));

disp('Running Demo1...');
Demo1;

disp('Running Demo2...');
Demo2;
close(fig4);

disp('Running Demo3...');
Demo3;

disp('All demos finished.');
