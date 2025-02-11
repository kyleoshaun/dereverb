close all
clear
clc

delete(gcp('nocreate'))

M = 6;
enable_source_whitening = true;

run_eval(M, enable_source_whitening)


