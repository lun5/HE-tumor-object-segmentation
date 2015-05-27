
addpath('/home/chakra/Bio/MEX-cca-yalla');

clear all
close all


load yaleInteraction

A = interactMat; 
% undirected graph loaded as upper triangular matrix
A = A + A';

% call to connected components on a sparse array  
[lbl z] = cca(A,0.1,1);

tt = find(lbl==1); 
B = A(tt,tt); 

figure; spy(B);

Kb = diag(sum(B,2)) - B;
[Uc,Sc] = eigs(sparse(Kb),50,'SA');
