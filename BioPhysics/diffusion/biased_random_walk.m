function [Xstore, Ystore]=biased_random_walk(I,N)
%function illustrates the biasd random work of N E coli within a solution
%of a log10(1/sqrt(x^2+y^2)) food density

%Input:     number of iterrations I
%           number of E. Coli N 



%all E Coli start randomly
xstart=-10;%+rand(1,N)*20;
ystart=-10;%+rand(1,N)*20;

xold=xstart;
yold=ystart;

%storing the path of all E Coli
Xstore=zeros(N,I);
Ystore=zeros(N,I);

Xstore(:,1)=xstart;
Ystore(:,1)=ystart;

for i=1:1:I
    
   
       
    %E coli has a 4sec memory where it measures the gradient
    if i/5-round(i/5)==0
        
        %location 4sec ago
        past_locationx=Xstore(:,i-4)';
        past_locationy=Ystore(:,i-4)';
        
        
        %measure gradient
        intermedx=0.5*(xnew+past_locationx);
        intermedy=0.5*(ynew+past_locationy);
        
        %this function determines the food density
        z1x=log10(1./sqrt(past_locationx.^2 + intermedy.^2)); 
        z2x=log10(1./sqrt(xnew.^2 + intermedy.^2)); 
        
        z1y=log10(1./sqrt(past_locationy.^2 + intermedx.^2));
        z2y=log10(1./sqrt(ynew.^2 + intermedx.^2)); 
        
        diffZx=(z2x-z1x);
        diffZy=(z2y-z1y);
        
        gradx=diffZx./(xnew-past_locationx);
        grady=diffZy./(ynew-past_locationy);
                     
        
        %Step in direction of the gradient. 
        
        sign_diffZxneg=find(gradx<0);
        sign_diffZyneg=find(grady<0);
                
        signx=ones(1,N);
        signx(sign_diffZxneg)=-1;      
       
        signy=ones(1,N);
        signy(sign_diffZyneg)=-1;     

        
        %If gradient is steep: step is
        %small in order to not miss the max. If gradient is flat: step is
        %large (but max ~five steps), since E Coli his hungry and likes 
        %to find food.
        xnew=xold+0.1*rand(1,N).*signx./(abs(gradx)+0.2);
        ynew=yold+0.1*rand(1,N).*signy./(abs(grady)+0.2);
        
        
        %no adaptive step size
        %xnew=xold+0.5*rand(1,N).*signx;
        %ynew=yold+0.5*rand(1,N).*signy;
                    
        
        
       
        
    else
        
        %pure random walk
    
        %random direction
        dix=rand(1,N);
        fx05=find(dix<0.5);
        dix(:)=1;
        dix(fx05)=-1;
    
        diy=rand(1,N);
        fy05=find(diy<0.5);
        diy(:)=1;
        diy(fy05)=-1;
    
    
        %do random step
    
        xnew=xold+0.1*rand(1,N).*dix;
        ynew=yold+0.1*rand(1,N).*diy;
       
        
    end
    
    %saving and resetting parameters
    
    Xstore(:,i)=xnew;
    Ystore(:,i)=ynew;
    
    xold=xnew;
    yold=ynew;
    
    
end







%for plotting
[X,Y]=meshgrid([-15:0.1:15], [-15:0.1:15]);
Z=log10(1./sqrt(X.^2 + Y.^2));
figure
contourf(X,Y,Z);
colormap gray

%plot individual paths of E Coli
   
 for n=1:1:N
         
     hold on     
     plot(Xstore(n,:),Ystore(n,:),'w-o')
     
 end
 %mark end and starting points
    hold on
    plot(Xstore(:,1),Ystore(:,1),'b+') 
    hold on
    plot(Xstore(:,end),Ystore(:,end),'r+')
    xlim([min(min(Xstore)) max(max(Xstore))])
    ylim([min(min(Ystore)) max(max(Ystore))])

if N>10
    
    %plotting histogram of effektive path length to the food source 
    diff=sqrt( Xstore(:,end).^2 + Ystore(:,end).^2);
    figure
    hist(diff,round(N/10));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','k')
    xlabel('distance from source')
    ylabel('#')

end














