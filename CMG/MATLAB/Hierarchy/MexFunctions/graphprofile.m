% function [B,D,M] = graphprofile(A,P);
%
% Input: adjaceny matrix A
%      : symmetric perturbation matrix P (optional)
%
% Output: D = sum(A)
%         M = max(A.*P)
%         B: uint32 vector representing unimodal tree 
%            such that M(i)=(A.*P)(B(i))

error ('graphprofile mexFunction not found') ;
