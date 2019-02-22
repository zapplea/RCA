% SYNTAX:
% m = RCA( 'arrayFilename', 'TFTargerList', 'zMode', nmfFlag, keepOriTFflag);
%
% INPUTS:
% arrayFilename - filename of array data matrix 
% TFTargerList - filename of targer gene list of TFs or miRNAs 
% zMode      - method to construct Z matrix. 'median', 'kmean', 'pls', and 'pca'
% nmfFlag   - binary flag. if set as '1' then Y is nonnegative 
%keepOriTFflag  - binary flag. if set as '1' keep original TF target


clear all;
close all;

m = RCA(‘input_array.txt', 'output1.txt', 'median',0, 1);
