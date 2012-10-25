close all; clear all;

addpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/

load sim_8ch_data

disp('data loaded');

inputImage = anatomy_orig;

numberOfSpokes = 64;

% generating the data vector

imageSize = size(inputImage);

theta = 0:numberOfSpokes-1;
theta = theta .* (180/numberOfSpokes);
dataMatrix = fftshift(fft2(fftshift(inputImage)));
dataMatrix = radon(dataMatrix, theta);
imageMatrix = iradon(dataMatrix, theta, 'linear', 'Ram-Lak', 1, imageSize(1));
imageMatrix = ifftshift(ifft2(ifftshift(imageMatrix)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.TVWeight = 0.77; 	% Weight for TV penalty
param.FOVWeight = 1;

% scale data
im_dc = imageMatrix;
dataMatrix = dataMatrix/max(abs(dataMatrix(:)));

im_dc = im_dc/max(abs(im_dc(:)));

res = im_dc;  %Initial degraded image supplied to fnlcg function
figure(300), imshow(abs(res), []);

% do iterations
tic
for n=1:5
	res = fnlCg(res,numberOfSpokes,dataMatrix, param);  %initialize fnlcg
	im_res = res;
	figure(100), imshow(abs(im_res),[]), drawnow;
end
toc


rmpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/
