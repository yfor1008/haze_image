close all
clear
clc

[file, path] = uigetfile('*.bmp; *.jpg; *.png', 'image ...');
im = imread([path file]);

% dehaze_im = dehaze(im, 7, 0.001, 0.95, 220, 0.1);
% dehaze_im_1 = imreducehaze(im, 'Method','approx','ContrastEnhancement','boost');
dehaze_im_1 = dehaze_approxDCP(im);
dehaze_im = dehaze(im);

figure, imshow(im)
figure, imshow(dehaze_im)
figure, imshow(dehaze_im_1)