function []=random_walk(I,xstart,ystart,zstart)

tic

x=xstart;
y=ystart;
z=zstart;
M=zeros(I+1,4);
M(1,:)=[x y z 0];
for i=1:1:I

    
    numb=rand;
    if numb>0.5
        sign_x=1;
    else
        sign_x=-1;
    end
    
    fx=sign_x*rand;
    
    numb=rand;
    if numb>0.5
        sign_y=1;
    else
        sign_y=-1;
    end
    
    fy=sign_y*rand;
    
    
    numb=rand;
    if numb>0.5
        sign_z=1;
    else
        sign_z=-1;
    end
    
    fz=sign_z*rand;

    x2=fx+x;
    y2=fy+y;
    z2=fz+z;
    
    %calculating current absolute distance from start point
    d=((x2-xstart)^2 + (y2-ystart)^2 + (z2-zstart)^2)^0.5;

    %saving values
    M(i+1,:)=[x2 y2 z2 d];
    
    %resetting values
    x=x2;
    y=y2;
    z=z2;
    
end

toc
plot3(M(:,1),M(:,2),M(:,3),'k-')
xlabel('x')
ylabel('y')
ylabel('z')
box on

figure
plot([1:1:I+1],M(:,4),'k-')
xlabel('iteration')
ylabel('distance from start point')




