% reads the input images...
% calls the motionEstNTSS function for further processing..

clear all,close all;
Image1=imread('18.tif');
Image2=imread('19.tif');
figure,imshow(Image1);
figure,imshow(Image2);

Image1gray=rgb2gray(Image1);
Image2gray=rgb2gray(Image2);

[motionVect, NTSScomputations,I] = motionEstNTSS(Image2gray, Image1gray, 8, 7);
