function [img1, img2] = readFramesGray(v1,v2)

img1 = readFrame(v1);
img2 = readFrame(v2);
img1 = rgb2gray(img1);
img2 = rgb2gray(img2);

