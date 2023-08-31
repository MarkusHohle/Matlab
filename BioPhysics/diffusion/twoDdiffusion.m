function [U]=twoDdiffusion() 


%set diffusion rate--------------------------------------------------------
Du=0.01;
%note: Du has to be small compared to one (=time step); otherwise the
%gaussian will decay into a ring (with central void) and then bounces back
%etc since the diffusion time within a numerical cell would be much less
%than a time step!

%define grid for XYT domains-----------------------------------------------
Lx=200;
Ly=200;
Lt=2000;
[X,Y] = meshgrid(1:1:Lx,1:1:Lx);                                
U=zeros(Lx,Ly,Lt);

%set initial conditions IC (gaussian)--------------------------------------
U(:,:,1)=round(2000*exp(-((X-round(Lx/2)).^2 + (Y-round(Ly/2)).^2)./round(Lx/50)));

%show IC ------------------------------------------------------------------

surf(U(:,:,1),'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong');
colormap(bone)
colorbar('location','southoutside')
zlim([0 2000]);
grid off
box on
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
        end
    end
    %plot every Xth time step
    if k-100*round(k/100)==0
        surf(U(:,:,k),'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong');
        colormap(bone)
        colorbar('location','southoutside')
        zlim([0 2000]);
        grid off
        box on
        pause
        close
    end
end


surf(U(:,:,k),'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','phong');
colorbar('location','southoutside')
colormap(bone)
grid off
box on
zlim([0 2000]);




















