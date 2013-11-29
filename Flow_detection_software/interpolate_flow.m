function interpolate_flow(path,mue2pix_ratio,fps,xcorr_thresh,k_size, k_sigma,k_size_temp,k_sigma_temp,max_flow_vel,flow_field_arrow_distance)
%Here we will load the retrograde flow data, and then interpolate it
%according to the spatial and temporal kernel size. The result will be
%stored in a subfolder. There will be a filtering: If the direction of a
%flow changes apruptly in time of space, it will be ignored.
filter_versioninfo='advanced_display_and_filter_v1_5';
%mue2pix_ratio=0.05
%fps=1/10;
%max_flow_vel=4; %this is in mue/min
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define kernel size. All kernel values have to be even numbers!!!!!!
if nargin < 4, xcorr_thresh=0.6; end                                                 %threshold value for cross-correlation
if nargin < 5, k_size=2*ceil(ceil(2*5/mue2pix_ratio)/2), end       %spatial gaussian filter kernel, the radius will be 5µm 
if nargin < 6, k_sigma=2*ceil(ceil(1/mue2pix_ratio)/2), end       %sigma of the spatial gaussian kernel should be 1µm
if nargin < 7, k_size_temp=4; end                                                    %time gaussian filter kernel, +- 2 frames
if nargin < 8, k_sigma_temp=2; end                                                %sigma of time gaussian
if nargin < 9, max_flow_vel=10; end                                                   %maximal flow velocity
if nargin < 10, flow_field_arrow_distance=25; end                               %distance between arrows in the interpolated images                    %maximal flow velocity


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Here we create the folder where we will later store the interpolation
%images and matlab files
mkdir([path,filesep,'filtered_overlay_im']);
mkdir([path,filesep,'filtered_overlay_im',filesep,'plus_actin']);
mkdir([path,filesep,'filtered_overlay_data']);


%first we will read the file list:
files=dir([path,filesep,'data',filesep,'*mat'])
for i=1:length(files)
    load([path,filesep,'data',filesep,files(i).name]);
    flow_data(i)=flow_struct;
    clear list
    %We directly switch from pixel per image to the mue/min
    list(:,1)=flow_struct.x_pos;
    list(:,2)=flow_struct.y_pos;
    list(:,3)=flow_struct.x_vec*mue2pix_ratio*fps*60;
    list(:,4)=flow_struct.y_vec*mue2pix_ratio*fps*60;
    list(:,5)=flow_struct.c_val;
    list_cell{i}=list;
end

%get the max and min x and y values of the rf vector startpoints for the full image series (means
%you get the values closest to the image boarder in all the time series)
t=max(size(list_cell));
for j=1:t
    list=list_cell{j};
    y(j)=max(list(:,2));
    x(j)=max(list(:,1));
    y_min(j)=min(list(:,2));
    x_min(j)=min(list(:,1));
end

%to prevent problems with rf values too close to the upper or left image boarder (so that the filter kernel would
%reach outside the image and negative array indices are not defined, in the positive direction there is no problem, the code
%just fills up the missing values with zeros),
%check if the min x,y values of the retroflow points in list_cell are within kernel/2.
%if not, we have to shift the startpoints of the flow arrows, and all the images in the series.
%all the shift info is stored in the x_shift, and y_shift variable.
if (min(x_min)<=k_size/2), x_shift=k_size/2-min(x_min)+1; else x_shift=0; end
if (min(y_min)<=k_size/2), y_shift=k_size/2-min(y_min)+1; else y_shift=0; end

%Save all parameters including x_shift and y_shift in a parameter file:
save([path,filesep,'analysis_parameters.mat'],'filter_versioninfo','path','xcorr_thresh','k_size','k_sigma','k_size_temp','k_sigma_temp','mue2pix_ratio','x_shift','y_shift');
    
%now shift the retroflow values for the full image series
for j=1:t
    list=list_cell{j};
    list(:,2)=list(:,2)+y_shift;
    list(:,1)=list(:,1)+x_shift;
    list_cell{j}=list;    
end

%define the k_size_temp+1 (e.g 5) layer thick convolution array block that is moved through the full
%image series in the analysis.
x=max(x);
y=max(y);
im_x(1:y+k_size+1+y_shift,2:x+k_size+1+x_shift,1:k_size_temp+1)=0;
im_devider=im_x;
im_y=im_x;

%generate the convolution kernel, gaussian kernel, circular, sigma of k_size/6
[k_x,k_y]=meshgrid(-k_size/2:+k_size/2,-k_size/2:+k_size/2);
kernel_2D=exp(-(k_x.^2+k_y.^2)/(2*(k_sigma)^2));
for i=1:k_size_temp/2+1
    k_size_temp-i+2;
    kernel_3D(:,:,i)=kernel_2D*exp(-(k_size_temp/2+1-i)^2/(k_sigma_temp)^2);
    kernel_3D(:,:,k_size_temp-i+2)=kernel_3D(:,:,i);    
end
min_val=max(min(kernel_3D(:,:,k_size_temp/2+1)));
kernel_3D(find(kernel_3D<min_val))=0;

startframe=1;
endframe=max(size(list_cell));
%endframe=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now that we have the unfiltered retroflow data, we can start the filtering 
%Do the 3D convolution
for j=startframe:endframe;
    list=list_cell{j};
    j;
    tic
    for i=1:max(size(list));
        x_c=list(i,1);
        y_c=list(i,2);
        dx=list(i,3);
        dy=list(i,4);
%        [th_d,r_d]=cart2pol(com_evolution(j,1)-x_c,com_evolution(j,2)-y_c);
        [th_f,r_f]=cart2pol(dx,dy);
%        th_diff=abs(th_d-th_f)*360/(2*pi);
%        if (th_diff>180), th_diff=abs(th_diff-360); end
        if(list(i,5)>xcorr_thresh) %&& th_diff<max_angle)
            x_c=list(i,1);
            y_c=list(i,2);
            lb_y=y_c-(k_size/2);
            ub_y=y_c+(k_size/2);
            lb_x=x_c-(k_size/2);
            ub_x=x_c+(k_size/2);
            %build the convolution matrix and the normalization matrix 
            im_x(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*list(i,3)*kernel_3D+im_x(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);
            im_y(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*list(i,4)*kernel_3D+im_y(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);
            im_devider(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:)=list(i,5)*kernel_3D+im_devider(y_c-(k_size/2):y_c+(k_size/2),x_c-(k_size/2):x_c+(k_size/2),:);          
        end
    end
    im_x_act=single(im_x(:,:,1)./im_devider(:,:,1));
    im_y_act=single(im_y(:,:,1)./im_devider(:,:,1));
    %create image and save final values (starts at frame k_size_temp/2 + 1 because of time filtering)
    if j>startframe-1+k_size_temp/2
        i=j-k_size_temp/2;
        %create the image, plus get the filtered final results
        create_retro_flow_image(path,flow_data(j).im,flow_data(j).edge_im,i,im_x_act,im_y_act,x_shift,y_shift,mue2pix_ratio,fps,max_flow_vel,flow_field_arrow_distance)
       % [final_res_filtered,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,i,im_x_act,im_y_act,final_res,x_shift,y_shift,mue2pix_ratio);
    end
    %delete first time layer of convolution array block and define next one
    im_x(:,:,1)=[];
    im_y(:,:,1)=[];
    im_devider(:,:,1)=[];
    im_x(:,:,k_size_temp+1)=0;
    im_y(:,:,k_size_temp+1)=0;
    im_devider(:,:,k_size_temp+1)=0;
    toc
    %create the last image save final values (has to be done differently because of time filtering)
    if j==max(size(list_cell));     %this means the last 2 frames run again but its easy this way
        for n=1:k_size_temp/2
            i=j-k_size_temp/2+n;
            im_x_act=single(im_x(:,:,n)./im_devider(:,:,n));
            im_y_act=single(im_y(:,:,n)./im_devider(:,:,n));
           % final_res=load([rf_path,filesep,'final_results', filesep,files_fin_res(i).name]);
            %create the image, plus get the filtered final results
            create_retro_flow_image(path,flow_data(j).im,flow_data(j).edge_im,i,im_x_act,im_y_act,x_shift,y_shift,mue2pix_ratio,fps,max_flow_vel,flow_field_arrow_distance)

%            [final_res_filtered,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,i,im_x_act,im_y_act,final_res,x_shift,y_shift,mue2pix_ratio);
        end
    end
 
end
close all;

    
