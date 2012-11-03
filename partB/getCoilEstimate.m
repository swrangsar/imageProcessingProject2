function coilEstimate = getCoilEstimate(inputImage, coilProfile, numberOfSpokes)


inverseCoilProfile = 1 ./ coilProfile;


% generating the data vector


theta = 0:numberOfSpokes-1;
theta = theta .* (180/numberOfSpokes);
dataMatrix = radon(inputImage, theta);
dataMatrix = fft(dataMatrix, [], 1);
imageMatrix = ifft(dataMatrix, [], 1);
imageMatrix = iradon(imageMatrix, theta, 'linear', 'Ram-Lak', 1, size(inputImage, 1));
imageMatrix = imageMatrix .* inverseCoilProfile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.TVWeight = 0.77;
param.FOVWeight = 1;
param.LaplacianWeight = 0.23;

param.CoilProfile = coilProfile;
param.InverseCoilProfile = inverseCoilProfile;

res = imageMatrix;
% figure, imshow(abs(res), []);
% title(['Coil data using ', num2str(numberOfSpokes), ' spokes']);


tic
figure;
for n=1:5
	[res, repetitionCounter] = fnlCgCoilEstimate(res,numberOfSpokes,dataMatrix, param);  %initialize fnlcg
	im_res = res;
	imshow(abs(im_res),[]), drawnow;
    if repetitionCounter > 5
        break;
    end;
end
title(['Coil estimate using ', num2str(numberOfSpokes), ' spokes']);
toc

coilEstimate = im_res;

end