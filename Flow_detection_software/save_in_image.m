function save_in_image(sig_x,sig_y,path,index)

    %This saves the arrays sig_x, sig_y in a png file, which is way better
    %compressed that the matlab data files!!

    abs_sig_x_min=abs(min(min(sig_x)));
    abs_sig_y_min=abs(min(min(sig_y)));
    sig_x_shift=sig_x+abs_sig_x_min;
    sig_y_shift=sig_y+abs_sig_y_min;
    sig_x_shift_max=max(max(sig_x_shift));
    sig_y_shift_max=max(max(sig_y_shift));
    sig_x_mat=mat2gray(sig_x_shift);
    sig_x_ind=gray2ind(sig_x_mat,65536);
    sig_y_mat=mat2gray(sig_y_shift);
    sig_y_ind=gray2ind(sig_y_mat,65536);
    imwrite(sig_x_ind,[path,filesep,'sig_x_ind_',num2str(index+1000),'_min_',num2str(abs_sig_x_min),'_shift_max_',num2str(sig_x_shift_max),'.png'],'png');
    filesep,imwrite(sig_y_ind,[path,filesep,'sig_y_ind_',num2str(index+1000),'_min_',num2str(abs_sig_y_min),'_shift_max_',num2str(sig_y_shift_max),'.png'],'png');
