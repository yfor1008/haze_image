function [B, T, atmLight] = dehaze_approxDCP(A)
% dehaze_approxDCP - ȥ��
%
% docs:
%   - matlabdd��imreducehaze.mʵ��
%

A = double(A) / 255; % [0, 1]
amount = 1;

% 1. Calculate dark channel image prior
patchSize = ceil(min(size(A,1), size(A,2)) / 400 * 15); % ��Сֵ�˲�����
minFiltStrel = strel('square', patchSize);

darkChannel = min(A,[],3); % �õ���ԭɫ
darkChannel = imerode(darkChannel, minFiltStrel); % ��ʴ����ȡ��Сֵ, ��Сֵ�˲�

% 2. Estimate atmospheric light
I = rgb2gray(A);
atmLight = estimateAtmosphericLight(A, I, darkChannel);

% 3. Estimate transmission t(x)
normI = A ./ atmLight;
normI = min(normI, [] , 3); % 
transmissionMap = 1 - imopen(normI, minFiltStrel); % ������(�ȸ�ʴ������), ȥ������, ��ʽ11

% 4. Use guided filtering to refine the transmission map
% Neighborhood size and degree of smoothing chosen
% empirically to approximate soft matting as best as possible.
epsilon = 1e-4;
filterRadius = ceil(min(size(A,1), size(A,2)) / 50);
nhoodSize = 2 * filterRadius + 1;
% Make sure that subsampleFactor is not too large
subsampleFactor = 4;
subsampleFactor = min(subsampleFactor, filterRadius);
transmissionMap = images.internal.algimguidedfilter(transmissionMap, A, [nhoodSize nhoodSize], epsilon, subsampleFactor); % ����ʵ�֣�����
transmissionMap = min(1, max(0, transmissionMap));
omega = 0.95;

% Thickness of haze in input image is second output of
% imreducehaze.Thickness Map does not depends on amount value.
T = 1 - transmissionMap;

% Omega value is set to 0.9, to leave some of haze in restored image for
% natural appearance of dehazed scene
transmissionMap = 1 - omega * (1 - transmissionMap); % ��ʽ12

% This lower bound preserves a small amount of haze in dense haze regions
t0 = 0.1;

% 5. Recover scene radiance
radianceMap = atmLight + (A - atmLight) ./ max(transmissionMap, t0); % ��ʽ22, t0����͸����̫С
radianceMap = min(1, max(0, radianceMap));

% New transmission map based on amount of haze to be removed
newTransmissionMap = min(1, transmissionMap + amount);

% Dehazed output image based on Amount, if Amount == 1,
% then B = radianceMap
B = radianceMap .* newTransmissionMap + atmLight .* (1-newTransmissionMap); %???

atmLight = double(reshape(atmLight, [1 3]));

% 6. Post-processing
B = boosting(B, amount, 1-T);

B = uint8(round(B * 255));
T = double(T);

end

function atmosphericLight = estimateAtmosphericLight(A, I, darkChannel)
% Atmospheric light estimation using 0.1% brightest pixels in darkchannel

% First, find the 0.1% brightest pixels in the dark channel.
% This ensures that we are selecting bright pixels in hazy regions.
p = 0.001; % 0.1 percent
[histDC, binCent] = imhist(darkChannel);
binWidth = mean(diff(binCent));
normCumulHist = cumsum(histDC)/(size(A,1)*size(A,2));
binIdx = find(normCumulHist >= 1-p);
darkChannelCutoff = binCent(binIdx(1)) - binWidth/2;

% Second, find the pixel with highest intensity in the
% region made of the 0.1% brightest dark channel pixels.
mask = darkChannel >= darkChannelCutoff;
grayVals = I(mask);
[y, x] = find(mask);
[~, maxIdx] = max(grayVals);

atmosphericLight = A(y(maxIdx(1)), x(maxIdx(1)), :);
atmosphericLight(atmosphericLight == 0) = eps(class(A));
end

function B = boosting(img, amount, transmissionMap)
% Boost as contrast enhancement technique
boostAmount = 0.1 * (1 - transmissionMap);
B = img .* (1 + (amount * boostAmount));
end

function enhanced = globalStretching(A)
% Global Stretching
chkCast = class(A);

% Gamma correction
gamma = 0.75;
A = A.^gamma;

% Normalization to contrast stretch to [0,1]
A = mat2gray(A);

% Find limits to stretch the image
clipLimit = stretchlim(A,[0.001, 0.999]);

% Adjust the cliplimits
alpha = 0.8;
clipLimit = clipLimit + alpha*(max(clipLimit, mean(clipLimit, 2)) - clipLimit);

% Adjust the image intensity values to new cliplimits
enhanced = imadjust(A, clipLimit);
enhanced = cast(enhanced, chkCast);

end