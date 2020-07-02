function [SegLabel,NcutDiscrete,NcutEigenvectors,W,imageEdges]= NcutImage(I,nbSegments,pts1,pts2,r,kappa)
%  [SegLabel,NcutDiscrete,NcutEigenvectors,NcutEigenvalues,W,imageEdges]= NcutImage(I);
%  Input: I = brightness image
%         nbSegments = number of segmentation desired
%  Output: SegLable = label map of the segmented image
%          NcutDiscrete = Discretized Ncut vectors
%  
% Timothee Cour, Stella Yu,  Jianbo Shi, 2004.


 
if nargin <2,
   nbSegments = 10;
end
    dataW.sampleRadius=r;
    dataW.sample_rate=1;
    dataW.edgeVariance = kappa;
[W,imageEdges] = ICgraph(I,dataW);

[nr,nc,nb] = size(I);

[B, c] = createB(pts1,pts2,nr,nc,W);
[NcutEigenvectors,NcutDiscrete] = ncut(W,nbSegments,B,c);
%% generate segmentation label map


SegLabel = zeros(nr,nc);
for j=1:size(NcutDiscrete,2),
    SegLabel = SegLabel + j*reshape(NcutDiscrete(:,j),nr,nc);
end