function [U1 U2]=diffusion_reaction() 
%reproduces fig 1c from Koch&Meinhard, Rev of Mod Phys, oct 1994
%pattern is stable after 500-600 time steps (press enter 10-12 times)

%set diffusion/depletion/production rates----------------------------------
Du1=0.005;
Du2=0.2;
rhou1=0.01;
rhou2=0.02;
ka=0.25;
%--------------------------------------------------------------------------

%define grid for XYT domains-----------------------------------------------
Lx=100;
Ly=100;
Lt=500;
%for u1
U1=zeros(Lx,Ly,Lt);
%for u2
U2=zeros(Lx,Ly,Lt);
%--------------------------------------------------------------------------

%set (random) initial conditions IC----------------------------------------
U1(:,:,1)=rand(Lx,Ly);
U2(:,:,1)=rand(Lx,Ly);
%plot ICs------------------------------------------------------------------
        figure
        %subplot(2,1,1)
        contourf(U1(:,:,1));
        colorbar('location','southoutside')
        colormap('gray')
        title('activator')
        %subplot(2,1,2)
        %surf(U1(:,:,1));
        %colorbar('location','southoutside')
        %colormap('gray')
        %title('activator')
        
        figure
        %subplot(2,1,1)
        contourf(U2(:,:,1));
        colorbar('location','southoutside')
        colormap('gray')
        title('inhibitor')
        %subplot(2,1,2)
        %surf(U2(:,:,1));
        %colorbar('location','southoutside')
        %colormap('gray')
        %title('inhibitor')
        
        fprintf('press ENTER to continue\n')
        
        pause
        close all
        
        fprintf('running....\n')
%--------------------------------------------------------------------------

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
                                                                                                                                                                                              %eps needed to avoid numerical artefacts 
            U1(i,j,k)=Du1*(U1(irun_up+1,j,k-1) + U1(irun_down-1,j,k-1) + U1(i,jrun_up+1,k-1) - 4*U1(i,j,k-1) + U1(i,jrun_down-1,k-1)) + U1(i,j,k-1) + rhou1*((U1(i,j,k-1)^2)/((1+ka*U1(i,j,k-1)^2)*U2(i,j,k-1)+eps) - U1(i,j,k-1));
            U2(i,j,k)=Du2*(U2(irun_up+1,j,k-1) + U2(irun_down-1,j,k-1) + U2(i,jrun_up+1,k-1) - 4*U2(i,j,k-1) + U2(i,jrun_down-1,k-1)) + U2(i,j,k-1) + rhou2*(U1(i,j,k-1)^2 - U2(i,j,k-1));
        end
    end
    
    %plot every fiftieth time step
    if k-50*round(k/50)==0
        figure
       % subplot(2,1,1)
        contourf(U1(:,:,k));
        colorbar('location','southoutside')
        colormap('gray')
        title('activator')
       % subplot(2,1,2)
       % surf(U1(:,:,k));
       % colorbar('location','southoutside')
       % colormap('gray')
       % title('activator')
        
        figure
       % subplot(2,1,1)
        contourf(U2(:,:,k));
        colorbar('location','southoutside')
        colormap('gray')
        title('inhibitor')
       % subplot(2,1,2)
       % surf(U2(:,:,k));
       % colorbar('location','southoutside')
       % colormap('gray')
       % title('inhibitor')
        
        fprintf('press ENTER to continue\n')
        
        pause
        close all
        
        fprintf('running....\n')
    end
end


%plot final state
        figure
      %  subplot(2,1,1)
        contourf(U1(:,:,k));
        colorbar('location','southoutside')
        colormap('gray')
        title('activator')
      %  subplot(2,1,2)
      %  surf(U1(:,:,k));
      %  colorbar('location','southoutside')
      %  colormap('gray')
      %  title('activator')
        
        figure
      %  subplot(2,1,1)
        contourf(U2(:,:,k));
        colorbar('location','southoutside')
        colormap('gray')
        title('inhibitor')
      %  subplot(2,1,2)
      %  surf(U2(:,:,k));
      %  colorbar('location','southoutside')
      %  colormap('gray')
      %  title('inhibitor')
















