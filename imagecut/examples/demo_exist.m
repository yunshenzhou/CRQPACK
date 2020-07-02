clear var;
close all;
addpath ../../src
addpath ../auxiliary
addpath ../data

%load image and labels from given examples, you may change the input to run
%different examples
[I,n1,n2,pts1,pts2,r,kappa] = imread_ncut('flower');


%compute the cut problem
tic;
disp('Start computing...');
nbSegments = 2;
[SegLabel,NcutDiscrete,NcutEigenvectors,W,imageEdges]...
     = NcutImage(I,nbSegments,pts1,pts2,r,kappa);
tt = toc;
disp(['The computation took ' num2str(tt) ' seconds on the '...
    num2str(size(I,1)) 'x' num2str(size(I,2)) ' image']);

%plot the labels
fig = figure(1);
clf;
imagesc(I);colormap(gray);axis off
hold on;
h = plot(pts1(:,1),pts1(:,2),'b+','MarkerSize',13);
set(h,'linewidth',5);
h = plot(pts2(:,1),pts2(:,2),'go','MarkerSize',13);
set(h,'linewidth',5);
hold off;
set(gca,'position',[0 0 1 1],'units','normalized');


%plot the result
fig = figure(2);
clf;
bw = edge(SegLabel,0.01);
J1 = showmask(I,imdilate(bw,ones(4,4))); imagesc(J1);axis off
hold on;
h = plot(pts1(:,1),pts1(:,2),'b+','MarkerSize',13);
set(h,'linewidth',5);
h = plot(pts2(:,1),pts2(:,2),'go','MarkerSize',13);
set(h,'linewidth',5);
hold off;
set(gca,'position',[0 0 1 1],'units','normalized');


%plot the eigenvector
fig = figure(3);
clf;
v = reshape(NcutEigenvectors,n1,n2);
imagesc(v);axis off;
set(gca,'position',[0 0 1 1],'units','normalized');

