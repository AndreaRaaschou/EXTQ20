function [t,y] = sim_equilibrium(t_end, y , p, dt )

k=p.g_d_vec;
A=p.A;
% disp_mat=p.disp_mat;
t=0; 




for i=1:t_end; 

    N=y(:,end);

    fN=(k*dt).*N+((A*dt)*N).*N; 
           
    y(:,end+1)=y(:,end)+fN;
    
    y(:,end)=max(y(:,end),0);
    
%     y(:,end)=disp_mat*y(:,end); %After the populations has reproduced and updated they disperse (Pontarp et al 2015)
        
    t(end+1)=t(end)+dt; 
    
end


