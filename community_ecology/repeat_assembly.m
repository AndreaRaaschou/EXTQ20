%repeat_assembly.m
%
%Lotka Volterra Assembly Competitive Communities
%
%This script calls the function
%
%final_num_species = assembly(c,s,num_invasions)
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


%Program control:
num_repetitions=100;
num_invasions=200;


%Vectors containing all variants values
%of c and amax we want to test

amax_vector=linspace(0,1,num_repetitions);
%c_vector=linspace(0,1,num_repetitions);
 

%Value of the parameter we do not test

c=0.5;
%amax=0.5;

%Main loop

for r=1:num_repetitions
   	amax=amax_vector(r);
    num_species_tot(r)=assembly(c,amax,num_invasions);
    
end

%End of Main loop

%Plot the number of species in final community against the varible we vary

figure(2)
plot(amax_vector,num_species_tot)
xlabel('Mean interaction strength')
%xlabel('Connectance')
ylabel('# Species in final community')



    