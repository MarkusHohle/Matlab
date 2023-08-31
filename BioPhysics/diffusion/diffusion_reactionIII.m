function [A S Y]=diffusion_reactionIII() 
%reproduces fig 9 from Koch&Meinhard, Rev of Mod Phys, oct 1994
%pattern is stable after 500-600 time steps (press enter 10-12 times)

%set diffusion/depletion/production rates----------------------------------
Da=0.01;
Ds=0.1;
rhoa=0.05;
rhos=0.0035;
rhoy=0.03;
mus=0.003;
muy=0.003;
sigs=0.0075;
sigy=0.00007;
ka=0.5;
ks=0.3;
ky=22;
%--------------------------------------------------------------------------

%define grid for XYT domains-----------------------------------------------
Lx=80;
Ly=80;
Lt=0.5e+4;
%for u1
A=zeros(Lx,Ly,Lt);
S=2.5*ones(Lx,Ly,Lt);
Y=zeros(Lx,Ly,Lt);
%--------------------------------------------------------------------------

%set initial conditions IC-------------------------------------------------
%uncommend the next three commended lines to add noise

%noise
A(:,:,1)=0.01*rand(Lx,Ly);
F=round(18*rand);
for f=1:1:F
    A(ceil(Lx*rand),ceil(Ly*rand),1)=2;
end
S(:,:,1)=S(:,:,1)+0.01*rand(Lx,Ly);
Y(:,:,1)=0.01*rand(Lx,Ly);

for k=2:1:Lt-1 %(we don't touch the ICs and Lt is calculated in the loop anyway)
    for j=1:1:Ly
        
         jrun_up=j;
         jrun_down=j;
         
         %cyclic BCs-------------------------------------------------------
            if j+1>Ly
                jrun_up=0;
            end
            if j-1==0
                jrun_down=Ly+1;
            end
          %----------------------------------------------------------------  
        
        for i=1:1:Lx 
            
            irun_up=i;
            irun_down=i;
            
            %cyclic BCs----------------------------------------------------
            if i+1>Lx
                irun_up=0;
            end
            if i-1==0
                irun_down=Lx+1;
            end
            %--------------------------------------------------------------
                                                                                                                                                                                               
            A(i,j,k)=Da*(A(irun_up+1,j,k-1) + A(irun_down-1,j,k-1) + A(i,jrun_up+1,k-1) - 4*A(i,j,k-1) + A(i,jrun_down-1,k-1)) + A(i,j,k-1) + rhoa*((A(i,j,k-1)^2)*S(i,j,k-1)/(1+ka*A(i,j,k-1)^2) - A(i,j,k-1));
            S(i,j,k)=Ds*(S(irun_up+1,j,k-1) + S(irun_down-1,j,k-1) + S(i,jrun_up+1,k-1) - 4*S(i,j,k-1) + S(i,jrun_down-1,k-1)) + S(i,j,k-1) + sigs/(1+ks*Y(i,j,k-1)) - rhos*S(i,j,k-1)*(A(i,j,k-1)^2)/(1+ka*A(i,j,k-1)^2) - mus*S(i,j,k-1);
            Y(i,j,k)=rhoy*(Y(i,j,k-1)^2)/(1+ky*Y(i,j,k-1)^2) + Y(i,j,k-1) - muy*Y(i,j,k-1) + sigy*A(i,j,k-1);  
        end
    end
    
    %plot every tenth time step
    if k-500*round(k/500)==0
        
        figure      
        contourf(Y(:,:,k));
        colormap(gray)
        colorbar('location','southoutside')
        fprintf('press ENTER to continue\n')
        pause
        close all
        fprintf('running....\n')
    end
end


%plot final state
        figure      
        contourf(Y(:,:,k));
        colormap(gray)
        colorbar('location','southoutside')















