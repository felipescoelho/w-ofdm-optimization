clc
clear
close all

mkdir('./logs')

diary ./logs/myFile.log

fprintf('ok\n')

diary off