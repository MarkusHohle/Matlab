function [U]=twoDdiffusion_smoluchouski(vx,vy) 
%input: components of flow vector
%eg: [U]=twoDdiffusion_smoluchouski(0.01,0.01);
%if  [U]=twoDdiffusion_smoluchouski(0.1,0.1); it shows a wavy structure:
%the flow is much faster then the diffusion

%Du and vx, vy have to be in the same order of magnitude to avoid waves

%set diffusion rate--------------------------------------------------------
Du=0.01;
%note: Du has to be small compared to one (=time step); otherwise the
%gaussian will decay into a ring (with central void) and then bounces back
%etc since the diffusion time within a numerical cell would be much less
%than a time step!

%define grid for XYT domains-----------------------------------------------
Lx=100;
Ly=100;
Lt=2000;
[X,Y] = meshgrid(1:1:Lx,1:1:Lx);                                
U=zeros(Lx,Ly,Lt);

%set initial conditions IC (gaussian)--------------------------------------
U(:,:,1)=round(2000*exp(-((X-round(Lx/2)).^2 + (Y-round(Ly/2)).^2)./round(Lx/50)));

%show IC ------------------------------------------------------------------

contourf(U(:,:,1));
colorbar('location','southoutside')
colormap(bone)
%zlim([0 2000]);
pause
close

for k=2:1:Lt
    for j=1:1:Ly
        
        jrun_up=j;
        jrun_down=j;
        %cyclic BCs--------------------------------------------------------
        if j+1>Ly
            jrun_up=0;
        end   
        if j-1==0
            jrun_down=Ly+1;           
        end
        %------------------------------------------------------------------ 
        
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
                                                                                                                                                                                         
                      %discrete laplacian 2D
            U(i,j,k)=Du*(U(irun_up+1,j,k-1) + U(irun_down-1,j,k-1) + U(i,jrun_up+1,k-1) - 4*U(i,j,k-1) + U(i,jrun_down-1,k-1))+ U(i,j,k-1);
                               %velocity times concentration gradient
            U(i,j,k)=U(i,j,k)-(vx/2)*(U(irun_up+1,j,k-1) - U(irun_down-1,j,k-1)) - (vy/2)*(U(i,jrun_up+1,k-1) - U(i,jrun_down-1,k-1));
        end
    end
    %plot every Xth time step
    if k-100*round(k/100)==0
        contourf(U(:,:,k));
        colorbar('location','southoutside')
        colormap(bone)
        %zlim([0 2000]);
        pause
        close
    end
end


surf(U(:,:,k));
colorbar('location','southoutside')
zlim([0 2000]);




















