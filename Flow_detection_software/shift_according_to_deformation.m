function I_out=shift_according_to_deformation(ux,uy,I_in)

%This moves the pixel values in I_in, according to ux,uy. If one pixel is reached multiple time, we average. 
%Intermediate pixels are filled with spline interpolation

ux=double(round(ux));
uy=double(round(uy));
I_in=double(I_in);

max_u=max(max(max(abs(ux))),max(max(abs(uy))));

[s_y,s_x]=size(ux);

sh_ux(1:s_y+2*max_u,1:s_x+2*max_u)=NaN;
sh_uy(1:s_y+2*max_u,1:s_x+2*max_u)=NaN;
sh_I(1:s_y+2*max_u,1:s_x+2*max_u)=NaN;
sh_ux(max_u+1:s_y+max_u,max_u+1:s_x+max_u)=ux;
sh_uy(max_u+1:s_y+max_u,max_u+1:s_x+max_u)=uy;
sh_I(max_u+1:s_y+max_u,max_u+1:s_x+max_u)=I_in;

ux=sh_ux;
uy=sh_uy;
I_in=sh_I;
I_in(isnan(I_in))=0;
[s_y,s_x]=size(ux);

if (size(ux)+size(uy)+size(I_in)~=3*size(ux) ) 
    'SIzes are not compatible'
    return
end

I_out(1:s_y,1:s_x,1)=0;
I_out(1:s_y,1:s_x,2)=NaN;
for i=1:s_x
    
    for j=1:s_y
        
        if(~isnan(ux(j,i)*uy(j,i)))
       %     I_out(j,i,1)=I_out(j,i,1)+I_in(j-uy(j,i),i-ux(j,i));
            I_out(j+uy(j,i),i+ux(j,i),1)=I_out(j+uy(j,i),i+ux(j,i))+I_in(j,i);
            if(isnan(I_out(j+uy(j,i),i+ux(j,i),2))) 
                    I_out(j+uy(j,i),i+ux(j,i),2)=1;
                
            else
                I_out(j+uy(j,i),i+ux(j,i),2)=I_out(j+uy(j,i),i+ux(j,i),2)+1;
            end
        end
    end
end

%Okay, to geth inner NaN's, we do a threshold, then close, then we try to
%detect the inner NaNs, and replace them by the mean of the "Non-NaNs" of
%the surrounding
SE = strel('square',3);
I_int=imerode(imdilate(I_out(:,:,1),SE),SE);
I_int2=I_out(:,:,2);
I_int(find(I_int==0))=-1;
I_int2(isnan(I_int2))=0;
I_int=I_int.*I_int2;


% Now go through again, and check if at the NaN posittions, you can get the
% value from the neighbours
kernel(1:3,1:3)=1;
kernel(2,2)=0;
kernel=kernel/8;
I_conv=conv2(I_out(:,:,1),kernel);
for i=1:s_x
    for j=1:s_y     
        if(I_int(j,i)==0)
             I_out(j,i,2)=1;
             I_out(j,i,1)=I_conv(j,i);
        end
    end
end

    
I_out=I_out(max_u+1:s_y-max_u,max_u+1:s_x-max_u,:);
I_out(:,:,3)=I_out(:,:,1)./I_out(:,:,2);

