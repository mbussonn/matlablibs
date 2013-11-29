function [final_res_out,min_pr,max_pr,max_r_flow]=advanced_create_il_v1_5(rf_path,index,x_flow,y_flow,final_res,x_shift,y_shift,mue2pix_ratio)
warning off MATLAB:divideByZero
hold off;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters (are saved in analysis_parameters.mat):
flow_field_arrow_distance=25;       %Flow field arrow distance
min_pr=0.1;                         %min value for polymerization grayscale
max_pr=6;                           %max value for polymerization grayscale
max_r_flow=5;                       %max retrograde flow amplitude for colorcoding
vel_edge_angles=[1:2:500];          %angles at which the edge velocity is drawn
max_edge_vel=10;                    %don't draw velocity bigger than max_edge_vel because these fluctuations are nonsense
stretch_vel_factor=10;              %multiplied with edge velocity to stretch the display arrows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first run: save the parameters defined above
if index==1
    save([rf_path,'\analysis_parameters.mat'],'flow_field_arrow_distance','min_pr','max_pr','max_r_flow','vel_edge_angles','max_edge_vel','stretch_vel_factor','-append');
end

%read the images 
files_im=dir([rf_path,'\movie\corr*']);
im1=imread([rf_path,'\movie\',files_im(index).name]);
[y_max,x_max]=size(im1);
%now take care of the x_shift, y_shift;
im_int(1:y_shift+y_max,x_shift+x_max)=0;
im_int(y_shift+1:y_shift+y_max,x_shift+1:x_shift+x_max)=im1;
im1=im_int;

%This is already prepared for Dannys Labview version, with the binary shapes stored in rf_path/binary_shape.
%If the folder exists, read from there, otherwise do it on your own
if isdir([rf_path,'\binary_shape']) 
    files_bin=dir([rf_path,'\binary_shape\corr*']);
    r_im_erode=double(imread([rf_path,'\binary_shape\',files_bin(index).name]));
    r_im_erode(find(r_im_erode==255))=1;
    r_im_erode_int(1:y_shift+y_max,x_shift+x_max)=0;
    r_im_erode_int(y_shift+1:y_shift+y_max,x_shift+1:x_shift+x_max)=r_im_erode;
    r_im_erode=r_im_erode_int;
else
    r_im_erode=better_edge(im1);    %edge detection function written by Danny.
end

[y_max,x_max]=size(im1);
max_size=max(size(im1));

x_flow=double(x_flow);
y_flow=double(y_flow);
[y_max_field,x_max_field]=size(x_flow);

%Make sure im1, and the flow data have the same size, if not reshape the flow data accordingly
if (~(size(x_flow)==size(im1)))
    intermediate(1:y_max,1:x_max)=NaN;
    intermediate(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]))=x_flow(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]));
    x_flow=intermediate;
    intermediate(1:y_max,1:x_max)=NaN;
    intermediate(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]))=y_flow(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]));
    y_flow=intermediate;
end
clear('intermediate');

%invert_im1
im1_inv=double(im1);
im1_inv=255-im1_inv;
im1_inv=uint8(im1_inv);

%first get the COM, and the edge values for the images, this is important to draw the edge velocity arrows,
%and the polymerization speed as grayscale value
se = strel('disk',1);
r_imlog=r_im_erode-imerode(r_im_erode,se);
[y_edge,x_edge]=find(r_imlog);
[r_y,r_x]=find(r_im_erode);
r_x_com=mean(r_x);
r_y_com=mean(r_y);
x_com=r_x_com;
y_com=r_y_com;
[th_edge,rad_edge]=cart2pol(x_edge-x_com,y_edge-y_com);
%use the shape, to cut of whatever of the retro flow maps is outside the real shape
x_flow=r_im_erode.*x_flow;
y_flow=r_im_erode.*y_flow;
x_flow(isnan(x_flow))=0;
y_flow(isnan(y_flow))=0;
[th_flow,r_flow]=cart2pol(x_flow,y_flow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DANNY: HERE I NEED TO SAVE THE INTERPOLATED FLOW AS PNG IMAGES TO SAVE DISK SPACE!!
%the available intensity range (65536) gives absolut precision to the 3rd after comma value (which is by far good enough)
dir_flow_path=[rf_path,'\dir_flow_analysis'];
%save([dir_flow_path,'\x_flow_array_',num2str(index+1000),'_dir_flow_analyis.txt'],'x_flow','-ASCII'); %old txt save. needs 1GB space
%save([dir_flow_path,'\y_flow_array_',num2str(index+1000),'_dir_flow_analyis.txt'],'y_flow','-ASCII'); %old txt save. needs 1GB space
abs_x_flow_min=abs(min(min(x_flow)));
abs_y_flow_min=abs(min(min(y_flow)));
x_flow_shift=x_flow+abs_x_flow_min;
y_flow_shift=y_flow+abs_y_flow_min;
x_flow_shift_max=max(max(x_flow_shift));
y_flow_shift_max=max(max(y_flow_shift));
x_flow_mat=mat2gray(x_flow_shift);
x_flow_ind=gray2ind(x_flow_mat,65536);
y_flow_mat=mat2gray(y_flow_shift);
y_flow_ind=gray2ind(y_flow_mat,65536);
%Reverse the x_shift and y_shift that was introduced in the "advanced_display_and_filter" routine
y_flow_backshift=y_flow_ind(y_shift+1:y_max,x_shift+1:x_max);
x_flow_backshift=x_flow_ind(y_shift+1:y_max,x_shift+1:x_max);
imwrite(x_flow_backshift,[dir_flow_path,'\x_flow_ind_',num2str(index+1000),'_min_',num2str(abs_x_flow_min),'_shift_max_',num2str(x_flow_shift_max),'.png'],'png'); 
imwrite(y_flow_backshift,[dir_flow_path,'\y_flow_ind_',num2str(index+1000),'_min_',num2str(abs_y_flow_min),'_shift_max_',num2str(y_flow_shift_max),'.png'],'png'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%now get the edge positions, r_vals, and th_vals in all the 500 directions
unsorted_edge_array(:,1)=x_edge;
unsorted_edge_array(:,2)=y_edge;
unsorted_edge_array(:,3)=round(th_edge/(2*pi)*500);
unsorted_edge_array(:,10)=th_edge/(2*pi)*500;
%to correct for the fact that the cart2pol makes negative angles, we have to find the negative ones and add 500
unsorted_edge_array(find(unsorted_edge_array(:,3)<=0),3)=unsorted_edge_array(find(unsorted_edge_array(:,3)<=0),3)+500;
unsorted_edge_array(find(unsorted_edge_array(:,10)<=0),10)=unsorted_edge_array(find(unsorted_edge_array(:,10)<=0),10)+500;
unsorted_edge_array(:,4)=rad_edge;

for i=1:max(size(y_edge)); 
    if unsorted_edge_array(i,3)==501
       unsorted_edge_array(i,3)=500;
    end
    unsorted_edge_array(i,5)=r_flow(y_edge(i),x_edge(i));
    unsorted_edge_array(i,6)=th_flow(y_edge(i),x_edge(i));
end

sorted_edge_array=sortrows(unsorted_edge_array,10);
x_edge_sorted_raw=sorted_edge_array(:,1);
y_edge_sorted_raw=sorted_edge_array(:,2);
clear('sorted_edge_array')
unsorted_edge_array(any(isnan(unsorted_edge_array)'),:) = [];

%The sorted_edge_array has is a 500x4 array:
%element(:,1)=angle with respect to COM, element(:,2)=distance to COM, element(:,3)=retro_flow_ampl, element(:,4)=retro_flow_angle
for i=1:500
    sorted_edge_array(1,i)=i;
    sorted_edge_array(2,i)=mean(unsorted_edge_array(find(unsorted_edge_array(:,3)==i),4));
    sorted_edge_array(3,i)=mean(unsorted_edge_array(find(unsorted_edge_array(:,3)==i),5));
    sorted_edge_array(4,i)=mean(unsorted_edge_array(find(unsorted_edge_array(:,3)==i),6));
end

%Now we need to take care of the directions that were not represented in the edge_image, or that had
%for some reasons no representatives
if (sum(isnan(sorted_edge_array(2,:)))>0)
    %there were detection problems, these are easily fixable for the distance, and retroflow amplitude.
    %However, for the angle we might run into inconsistensies
    sorted_edge_array(:,find(isnan(sorted_edge_array(2,:))))=[];
    sorted_edge_array_interp(1,:)=[1:500];
    sorted_edge_array_interp(2,:)=interp1(sorted_edge_array(1,:),sorted_edge_array(2,:),[1:500]);
    sorted_edge_array_interp(3,:)=interp1(sorted_edge_array(1,:),sorted_edge_array(3,:),[1:500]);
    sorted_edge_array_interp(4,:)=interp1(sorted_edge_array(1,:),sorted_edge_array(4,:),[1:500]);
    sorted_edge_array=sorted_edge_array_interp;
end

%So, sorted_edge_array represents now the com-edge distance, the retroflow, and the direction of retro flow in each direction.
%Now we just need to combine this information with the angular velocities that we got from the final results array

%In the current labview program the velocity is defined positive for outward motion.
%To be in line with the old analysis, we need to check, if the files were analyzed with the old program
if exist([rf_path,'\Readme.txt'])==0    %if Readme.txt not there it's old analyis and needs multiply by -1.
    velocity(1,:)=-1*final_res(:,2)';   %old definition: values have to be multiplied by -1
else
    velocity(1,:)=final_res(:,2)';      %standard now
end
velocity(2,:)=final_res(:,4)';

%reassemble the filtered final_res, out of the data we got so far.
%Additionally to the 4 fields: retro_flow ampl, edge_vel, edge_pol, relyability factor.
%we add retro_flow direction, corrected (for retroflow direction) polymerization speed, com_edge distance.
final_res_out(:,1)=sorted_edge_array(3,:)';                     %Retroflow amplitude
final_res_out(:,2)=velocity(1,:)';                              %Edge velocity   
final_res_out(:,3)=final_res_out(:,1)+final_res_out(:,2);       %Edge Polymerization = Retroflow + Edge Velocity
final_res_out(:,4)=velocity(2,:)';                              %Reliability factor
final_res_out(:,5)=sorted_edge_array(4,:)';                     %Retroflow direction
final_res_out(:,6)=final_res_out(:,2)+final_res_out(:,1).*cos((sorted_edge_array(1,:)+250)'/500*2*pi-sorted_edge_array(4,:)');  %Polymerization Speed corrected for retroflow vector angle (the +250/500*2*pi is to convert angle to standard coordinate system)

%since we might have a lot of badly detected edge velocity values, also interpolate the stupid ones!!
speed=[1:500]';
speed(:,2)=velocity(1,:)';
speed(:,3)=velocity(2,:)';
speed(find(~speed(:,3)),:)=[];
interp_speed(:,1)=[1:500]';
interp_speed(:,2)=interp1(speed(:,1),speed(:,2),[1:500]');

%Now we assemble everything to get the display_edge_array, which has the form:
%1. angle(1:500),2. retro_flow_ampl, 3. interpol_speed, 4. Recalc polymerization rate, 5. edge_distance, 6. retro_flow_angle,
%7. Recalc_polymerization rate using retro_flow_angle 8. x-component of according edge point, 9. y-component of according edge point
display_edge_array=[1:500]';
display_edge_array(:,2)=final_res_out(:,1);
display_edge_array(:,3)=interp_speed(:,2);

%smooth the velocity for the display_edge_array
edge_vel=display_edge_array(:,3);
edge_vel(501:1000)=display_edge_array(:,3);
edge_vel(1001:1500)=display_edge_array(:,3);
edge_vel=smooth(edge_vel);
display_edge_array(:,3)=edge_vel(501:1000);

%Stretch the velocity amplitude for the later display by the strech factor
display_edge_array(:,3)=stretch_vel_factor*display_edge_array(:,3);
clear('edge_vel');

display_edge_array(:,4)=display_edge_array(:,3)/stretch_vel_factor+display_edge_array(:,2);
display_edge_array(:,5)=sorted_edge_array(2,:)';
display_edge_array(:,6)=sorted_edge_array(4,:)';
display_edge_array(:,7)=display_edge_array(:,3)/stretch_vel_factor+display_edge_array(:,2).*cos((sorted_edge_array(1,:)+250)'/500*2*pi-sorted_edge_array(4,:)');
display_edge_array(:,8)=x_com+display_edge_array(:,5).*cos(sorted_edge_array(1,:)'/500*2*pi);
display_edge_array(:,9)=y_com+display_edge_array(:,5).*sin(sorted_edge_array(1,:)'/500*2*pi);

%Now we should have everything together to start the drawing the image
%First make the color coded growth retrograde flow
r_flow(isnan(r_flow))=0;
imshow(r_flow)
colormap('jet')
set(gca,'CLim',[0,6])
hold on

%now get and draw the flow field arrows
plot_data=generate_plot_normalized(x_flow,y_flow,flow_field_arrow_distance);
hold on
quiver(plot_data(:,1),plot_data(:,2),plot_data(:,3),plot_data(:,4),0.3,'y');
hold off

%next, redraw the edge, by giving the grayscale colorcoding for the polymerization amplitude

%now assign the pol val to the x,y edge coordinates by abusing the unfiltered_edge_array
%unsorted_edge_array(i,7):polymerization rate without using the retro angle
%unsorted_edge_array(i,8):polymerization rate with using the retro angle
for i=1:max(size(unsorted_edge_array)); 
    unsorted_edge_array(i,7)=display_edge_array(unsorted_edge_array(i,3),4);
    unsorted_edge_array(i,8)=display_edge_array(unsorted_edge_array(i,3),7);
end
sorted_edge_array=sortrows(unsorted_edge_array,10);

%Now interpolate the values for the edge, to create a smooth colorcoding in the grayscale
z=griddata(sorted_edge_array(:,1),sorted_edge_array(:,2),sorted_edge_array(:,8),x_edge_sorted_raw,y_edge_sorted_raw);
z=round((z-min_pr)*127/(max_pr-min_pr)+1);
z(z>=127)=127;
z(z<1)=1;
hold on
cmp=gray(128);
n=max(size(x_edge_sorted_raw));
for i=1:n
    if isnan(z(i))
        if i==1
            z(i)=z(500);
        else
            z(i)=z(i-1);
        end
    end
    %now check if the z(i) is still NaN, if so the polymerization plot breaks down
    %-> report this in the frame, and stop the polymerization plot
    if isnan(z(i))
        text(50,10,'Polymerization interpolation failed!!!','Color','r')
        break;
    end
    plot(x_edge_sorted_raw(i),y_edge_sorted_raw(i),'LineStyle','none','Marker','.','MarkerSize',5,'MarkerEdgeColor',cmp(z(i),:))
end

%Draw a 10µm scale bar
fill([x_max-5-10/mue2pix_ratio,x_max-5-10/mue2pix_ratio,x_max-5,x_max-5],[5,15,15,5],'w')
%text(x_max-5-10/mue2pix_ratio,20,'10 µm','Color','w')

%before drawing the edge velocity arrows save a version without the edge velocity arrows, to make a nice movie
out=[rf_path,'\filtered_overlay_im\without_edge_vel','\all_overlay',num2str(index+1000),'all_filt.png'];
print('-dpng','-r150',out);

%And last but not least, draw the edge velocity arrows. To get rid of strong ones, smooth the array
display_edge_array(find(abs(display_edge_array(:,3))>max_edge_vel*stretch_vel_factor),3)=NaN; %Here we look for edge velocites that are bigger than the max_edge_vel
v=vel_edge_angles;
quiver(display_edge_array(v,8),display_edge_array(v,9),display_edge_array(v,3).*cos((display_edge_array(v,1))*2*pi/500),display_edge_array(v,3).*sin((display_edge_array(v,1))*2*pi/500),0,'r');

%Draw a little 10µm/min arrow to see velocity scale
quiver(10,10,0,max_edge_vel*stretch_vel_factor,'r')
text(20,10,[num2str(max_edge_vel),'µm/min'])

%Save the images, the display_edge_array, and of course the final_results_out

out=[rf_path,'\filtered_overlay_im','\all_overlay',num2str(index+1000),'all_filt.png'];
print('-dpng','-r150',out);
%now correct for the stupid matlab exporting and resizeing
im_corr=imread(out);
im_bw=mean(im_corr,3);
x_mean=mean(im_bw,2);
y_mean=mean(im_bw,1);
x_im_f=find(~(x_mean==255));
y_im_f=find(~(y_mean==255));
im_corr_cut=im_corr(x_im_f(1):x_im_f(length(x_im_f)),y_im_f(1):y_im_f(length(y_im_f)),:);
[y_im_corr_cut,x_im_corr_cut,depth]=size(im_corr_cut);
im1_resized=imresize(im1_inv,[y_im_corr_cut,x_im_corr_cut]);
imwrite(im_corr_cut,out);
%finally get an overlay with the original image, to give the structural information!!
im1_inv_d(:,:,1)=double(im1_resized);im1_inv_d(:,:,2)=double(im1_resized);im1_inv_d(:,:,3)=double(im1_resized);
im_corr_cut_im1=uint8(ceil(double(im_corr_cut).*im1_inv_d./255));
imwrite(im_corr_cut_im1,[rf_path,'\filtered_overlay_im','\plus_actin\all_overlay',num2str(index+1000),'all_filt.png']);
save([rf_path,'\filtered_final_res','\filtered_',num2str(index+1000),'_final_res.txt'],'final_res_out','-ASCII'); 
save([rf_path,'\filtered_overlay_im\display_edge_array','\filtered_',num2str(index+1000),'display_edge_array.txt'],'display_edge_array','-ASCII');

toc


