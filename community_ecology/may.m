%may.m
%
%Stability and complexity for ecological communities
%
%This script calls the function
%
%   J = random_matrix(s,n,c)
%
%
%
%Excercise Theoretical Ecology
%Lund November 2005
%Jacob Johansson
%
%


%Clear the workspace and close figures
clear all
close all


%Parameters
%Interaction strength standard deviation
s=0.2;
%Connectance
c=0.5;

%Number of repetitions
num_repetitions=200;

%The different values for the variable we change
%is stored in a vector which is num_repetition long 

n_vector=1:num_repetitions;


for r=1:num_repetitions
    
    %Pic a value for
    %s=s_vector(r);
    n=n_vector(r);

    
    %Create random matrix
    J=random_matrix(s,n,c);
    
    %Find the dominant eigenvalue
    lambda_max=max(real(eig(J)));

    %Record the dominant eigenvalue for plotting later
    lambda_max_tot(r)=lambda_max;
    
    %Record the s *sqrt(nc) for plotting later
    snc(r)=s*sqrt(n*c);


end


%Plot the dominant eigenvalue against s*sqrt(nc)
plot(snc,lambda_max_tot,'k.')
xlabel('s\surdnc')
ylabel('dominant eigenvalue')