function [I,Inr,Inc,pts1,pts2,r,kappa] = imread_ncut(imgName)
%  I = imread_ncut(Image_file_name);
%
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.


%% read image 
I = imread(strcat('../data/', imgName, '.jpg'));
[Inr,Inc,nb] = size(I);

if (nb>1),
    I =double(rgb2gray(I));
else
    I = double(I);
end

%load points
load(strcat('../data/pts_', imgName, '.mat'));
pts1 = cell2mat(pts(1));
pts2 = cell2mat(pts(2));

%set parameters
switch imgName
    case 'crab'
        kappa = 0.1;
        r = 10;
    case 'camel'
        kappa = 0.08;
        r = 5;
    case 'face1'
        kappa = 0.1;
        r = 5;
    case 'face2'
        kappa = 0.1;
        r = 5;
    case 'daisy'
        kappa = 0.08;
        r = 5;
    case 'daisy2'
        kappa = 0.08;
        r = 5;
    case 'flower'
        kappa = 0.1;
        r = 5;
    otherwise
        kappa = 0.1;
        r = 5;
end
