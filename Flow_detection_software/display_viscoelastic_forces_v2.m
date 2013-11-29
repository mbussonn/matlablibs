function display_viscoelastic_forces_v2(path)

%This function just reads the forces stored in png images, converts them to
%polar coordinates, displays them, and finally stores them

Max_F=30;
flow_field_arrow_distance=20;


max_i=length(sig_x_im_name_in);

for index=1:max_i
    index
    xstring=sig_x_im_name_in(index).name;
    ystring=sig_y_im_name_in(index).name;
    sig_x_min=str2num(xstring(findstr('min_',xstring)+4:findstr('_shift',xstring)-1));
    sig_y_min=str2num(ystring(findstr('min_',ystring)+4:findstr('_shift',ystring)-1));
    sig_x_shift_max=str2num(xstring(findstr('max_',xstring)+4:findstr('.png',xstring)-1));
    sig_y_shift_max=str2num(ystring(findstr('max_',ystring)+4:findstr('.png',ystring)-1));

    sig_x_im_in=imread([path,'\internal_forces\',sig_x_im_name_in(index).name]);
    sig_y_im_in=imread([path,'\internal_forces\',sig_y_im_name_in(index).name]);

    sig_x_im_back=double(sig_x_im_in)./65536.*sig_x_shift_max;
    sig_y_im_back=double(sig_y_im_in)./65536.*sig_y_shift_max;
    sig_x=sig_x_im_back-sig_x_min;
    sig_y=sig_y_im_back-sig_y_min;
    
    sig_x(find(sig_x==sig_x(1,1)))=NaN;
    sig_y(find(sig_y==sig_y(1,1)))=NaN;
    
    [sig_th,sig_r]=cart2pol(sig_x,sig_y);
    
    %now start and make the plot!

    %First plot the radial force part
    h=figure;
    imshow(sig_r)
    colormap('jet')
    set(gca,'CLim',[0,Max_F])
    colorbar
    hold on

    %now get the flow field arrows
    plot_data=generate_plot_normalized(sig_x,sig_y,flow_field_arrow_distance);
    quiver(plot_data(:,1),plot_data(:,2),plot_data(:,3),plot_data(:,4),0.3,'y');
    print('-dpng','-r300',[path,'\internal_forces\','display_vis_',num2str(index+1000),'.png']);
    %saveas(h,[path,'\internal_forces\',num2str(index),'.fig']);
    hold off
    close
    max(max(sig_r))
    

end
