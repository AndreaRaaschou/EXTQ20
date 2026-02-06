%assembly.m
%
%Lotka Volterra Assembly Competitive Communities
%
%
%Excercise Theoretical Ecology
%Lund November 2005
%Jacob Johansson
%
%
%function final_num_species=assembly(c,s,num_invasions)
%
%Parameters in:
%
%Connectance:                   c
%Mean interaction strength:     s
%How many invasion attempts:    num_invasions
%
%
%
%This function calls the following functions:
%
%
%n_star         =   EqDensity_Func(A);
%[arow,acol]    =   NewSpecies_Func(no_spec,c,s);
%A              =   NewMatrix_Func(A,Nstar);
%                   InvGrowth_Func(arow,n_star)



function final_num_species=assembly(c,amax,num_invasions)


%Initial settings
A=-1;
num_spec=1;



%Main loop

for t=1:num_invasions

    %Calculate equilibrium densities
    n_star=EqDensity_Func(A);

    %Draw a new species from the pool
    [a_row,a_col]=NewSpecies_Func(num_spec,c,amax);

    %Test if the invading species can grow
    %initially in the new community
    if InvGrowth_Func(a_row,n_star)

        %Add the new species to the community
        A=[A;a_row];
        A=[A,a_col];
        num_spec=num_spec+1;


        %Calculate equilibrium densities including the invading species
        n_star=EqDensity_Func(A);

        %Remove negative populations until all positive
        while any(n_star<0)

            %Remove corresponding interactions from the matrix
            A=NewMatrix_Func(A,n_star);

            %Calculate the new species number
            num_spec=size(A,1);

            %Calculate equilibrium densities when
            %negative densities have been removed
            n_star=EqDensity_Func(A);


        end
    end
    
    %Save the number of species in a vector
    %in order to plot it later
    num_species_tot(t)=num_spec;

end



%End of Main loop



%Plot the number of species over time
figure(1)
plot(1:num_invasions,num_species_tot)
title('')
xlabel('# Invasion attempts')
ylabel('# Species')

%Give the function output
final_num_species=num_spec;


