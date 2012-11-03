close all; clear all;

addpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/

load sim_8ch_data;
disp('sim_8ch_data loaded successfully!');


inputImage = data;
coilProfile = b1;

numberOfSpokes = 32;
numberOfCoils = 8;

coilEstimate = cell(numberOfCoils, 1);

for k = 1:numberOfCoils
    coilEstimate{k} = getCoilEstimate(data(:, :, k), b1(:, :, k), numberOfSpokes);
%     figure(k*100), imshow(abs(coilEstimate{k}), []);
end

avgCoilEstimate = zeros(size(coilEstimate{1}));
for k = 1:numberOfCoils
    avgCoilEstimate = avgCoilEstimate + coilEstimate{k};
end
avgCoilEstimate = avgCoilEstimate ./ numberOfCoils;

% since we removed the phases we can now take just the real part of the avg
% estimate image

finalCoilEstimate = real(avgCoilEstimate);


% figure(900); clf; imshow(abs(avgCoilEstimate), []);

%% next we improve the average coil estimate

% generating the data vector

theta = 0:numberOfSpokes-1;
theta = theta .* (180/numberOfSpokes);
dataMatrix = radon(finalCoilEstimate, theta);
dataMatrix = fft(dataMatrix, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.TVWeight = 0.0001; 	% Weight for TV penalty
param.FOVWeight = 10;
param.POSWeight = 5;
param.LaplacianWeight = 0.23;



res = finalCoilEstimate;  %Initial degraded image supplied to fnlcg function
figure(300), imshow(abs(res), []);

% do iterations
tic
for n=1:5
	[res, repetitionCounter] = fnlCg(res,numberOfSpokes,dataMatrix, param);  %initialize fnlcg
	im_res = res;
	figure(100), imshow(abs(im_res),[]), drawnow;
    title(['Image estimate using ', num2str(numberOfSpokes), ' spokes from multiple coils']);

    if repetitionCounter > 7
        break;
    end;
end
toc

rmpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/