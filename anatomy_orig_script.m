close all; clear all;

addpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/

load sim_8ch_data

disp('data loaded');

inputImage = anatomy_orig;

numberOfSpokes = 8;

% generating the data vector


theta = 0:numberOfSpokes-1;
theta = theta .* (180/numberOfSpokes);
dataMatrix = radon(inputImage, theta);
dataMatrix = fft(dataMatrix, [], 1);
imageMatrix = ifft(dataMatrix, [], 1);
imageMatrix = iradon(imageMatrix, theta, 'linear', 'Ram-Lak', 1, size(inputImage, 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.TVWeight = 0.0001; 	% Weight for TV penalty
param.FOVWeight = 1;
param.POSWeight = 5;
param.LaplacianWeight = 0.23;



res = imageMatrix;  %Initial degraded image supplied to fnlcg function
figure(300), imshow(abs(res), []);
title(['Initial image estimate using ', num2str(numberOfSpokes), ' spokes']);

% do iterations
tic
for n=1:5
	[res, repetitionCounter] = fnlCg(res,numberOfSpokes,dataMatrix, param);  %initialize fnlcg
	im_res = res;
	figure(100), imshow(abs(im_res),[]), drawnow;
    title(['Image estimate using ', num2str(numberOfSpokes), ' spokes']);
    
    
    if repetitionCounter > 7
        break;
    end;
end
toc


rmpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/