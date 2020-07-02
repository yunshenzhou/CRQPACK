clear var;
close all
addpath ../../src
addpath ../auxiliary
addpath ../data

%% read image 
imgName = 'bird';
I = imread(strcat( imgName, '.jpg'));
I2 = imresize(I,1);
[n1,n2,nb] = size(I2);

if (nb>1)
    I  = double(rgb2gray(I));
else
    I = double(I);
end
%set the paremeters
r = 5;
kappa = 0.1;

%% label the points
figure(1);
clf; 
imagesc(I);
colormap(gray);
axis off;
set(gca,'position',[0 0 1 1],'units','normalized');
disp('Click the labels of the background and press Enter to finish');
[x, y] = ginput;
pts1 = [x,y];
hold on;
h = plot(pts1(:,1),pts1(:,2),'b+','MarkerSize',13);
set(h,'linewidth',5);
disp('Click the labels of the object and press Enter to finish');
[x, y] = ginput;
pts2 = [x,y];
h = plot(pts2(:,1),pts2(:,2),'go','MarkerSize',13);
set(h,'linewidth',5);
% save the labels
%save(strcat('../data/pts_',imgName,'.mat'),'pts1','pts2');

%%compute the cut
disp('Start computing...');
nbSegments = 2;
tic;
[SegLabel,NcutDiscrete,NcutEigenvectors,W,imageEdges]...
     = NcutImage(I,nbSegments,pts1,pts2,r,kappa);
tt = toc;
disp(['The computation took ' num2str(tt) ' seconds on the '...
    num2str(size(I,1)) 'x' num2str(size(I,2)) ' image']);

%% plot the results
figure(2);
clf;
bw = edge(SegLabel,0.01);
J1 = showmask(I,imdilate(bw,ones(1,1))); imagesc(J1);axis off
hold on;
h = plot(pts1(:,1),pts1(:,2),'b+','MarkerSize',13);
set(h,'linewidth',5);
h = plot(pts2(:,1),pts2(:,2),'go','MarkerSize',13);
set(h,'linewidth',5);
hold off;
set(gca,'position',[0 0 1 1],'units','normalized');


