function [B, c] = createB(pts1,pts2,nr,nc,W)
%
%   purpose:
%      create the matrix B for the constraints Bv=0
%   Input:
%      pts1: points selected and labeled for cluster 1
%      pts2: points selected and labeled for cluster 2
%      nr:   number of rows of image I
%      nc:   number of coumns of image I
%   Output:
%      B: the constraints matrix
%
n = nr*nc;
BP1 = round(pts1);
BP2 = round(pts2);
x1 = BP1(:,1);
y1 = BP1(:,2);
x2 = BP2(:,1);
y2 = BP2(:,2);
m1 = size(BP1,1);
m2 = size(BP2,1);
%%   homogenous
% B = zeros(m1+m2-1,n);
% idx1 = (x1(1)-1)*nr+y1(1);
% for i = 2:m1
%     idx = (x1(i)-1)*nr+y1(i);
%     B(i-1, idx1)=1;
%     B(i-1, idx)= -1;
% end
% for i = m1+1:m1+m2
%     idx = (x2(i-m1)-1)*nr+y2(i-m1);
%     B(i-1, idx1)=1;
%     B(i-1, idx) =1;
% end
% c = zeros(m1+m2-1,1);
%% imhomogenous
B = zeros(m1+m2,n);
c = zeros(m1+m2,1);
for i = 1:m1
    %j = (x1(i)-1)*nr+y1(i)
    j = (x1(i))+(y1(i)-1)*nc;
    j = (x1(i)-1)*nr+y1(i);
    B(i,j) = 1;
    c(i) = 1;
end
for i = 1:m2
    j = (x2(i))+(y2(i)-1)*nc;
    j = (x2(i)-1)*nr+y2(i);
    B(i+m1,j) = 1;
    c(i+m1) = -1;
end









