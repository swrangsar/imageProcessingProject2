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

figure(900); clf; imshow(abs(avgCoilEstimate), []);

rmpath /Users/swrangsarbasumatary/Desktop/imageProcessingProject2/