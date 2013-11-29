function Pol=get_viscoelastic_forces_v2(path)
%using the displacement fields, we calculate the stresses that are required
%to produce these deformtions. This is done by assuming a stress relaxation
%function that can be either modeled by a exponential decay, or be measured
%with an AFM. The viscoelastic responce is the nothing but the decay
%function multiplied with the Elasticity coefficient in the of the pure
%elastic case, multiplied this the strain rate, integrated over time!

%The mathematical model is described in the supplemental material of:
%Betz, "Growth cones as soft and weak force generators." PNAS 2011 Aug 16;108(33)

% Okay, we do this the folowing way:
% First read the first deformations, and calculate the according stresses (in Array S),
% then shift them, and multiply with exp(-delta_t/tau).
% After that, do the next images, and put in it in S(:,:,2), the sum along
% the 3rd dimension is the resulting stress!
% After that repeat the shift and the relaxation. Do this over and over
% until done. 

%First we need to get the path where the flow fields are stored.
%The program assume that in ech file there is one dataset, where the x and y flow 
%are stored under the name x_flow y_flow, each a 2 array that represents
%the flow value for each pixel in the original Image. Furthemore the frame
%per second is stored under fps, and the ratio between mue to pixel is
%stored under mue2pix_ratio
%If the path is not given the program will ask for it
if nargin < 1, 
    path=uigetdir('E:\','Give me the basepath with the Info file of the data to analyze.');
end
%path='E:\Science\data\brian_stemer\test_temp1\filtered_overlay_data'


Max_F=200;%This parameter is the maximal force for the colormap of the final images
flow_field_arrow_distance=20 %THsi parameter gives the distance between flow field arrows in the final image


%path='D:\Data_retr\GFP-Actin NG108 050205 Big Stationary [done,++]\GFP-Actin NG108 050205 Big Stationary [done,++]\standard\retrograde_flow_data';%If you change this, change the mue2pix ratio too


%Here we define the coefficients we got out of the AFM-measurements for the
%growth cone, using the dashpot - spring||dashpot model. These values are
%measured with YunBe. They correspond to: E=112, tau=2.4, eta=4.8
%These values need to be remeasured for each cell type used.
a1=905.7536;
a2=17.6628;
b1=8.5519

%load the first dataset
flow_files=dir([path,'\*.mat']);
load([path,'\',flow_files(3).name]);
%Then we try to extract some important parameters:
mue2pix=mue2pix_ratio;
delta_t=1/fps;


%According to the dashpot- spring||dashpot model, the prefactor before the
%exponential decay reads:
E=(a1*b1-a2)/b1^2;
%NOw we convert it to pixel units:
E=E*mue2pix^2;
%Poisson Ratio
g=0.47;


%relaxation time
tau=b1;
%since the measured deformation always takes place during the delta_t time
%interval, we model the strain rate to be constant: ux/delta_t
%Putting this together with the assumption of the exponential stress
%relaxation, we can actually calculate the int exp(t/tau)*ux/delta_t dt.
%By devinding through the delta_t, we already include time required for the
%strain in this coefficient. Thus We can just use the pure elastic
%calculation!
%This gives
int_decay=tau*(1-exp(-delta_t/tau));

%Now the effective E is nothing but the prefactor from the dashpot spring
%dashpot model, times this "int_decay" + the constant term in that results
%from the model a2/b1 timesw the delta_t (dont forget to convert the A2/b1
%to pixel units!!
E=E*int_decay+0.5*a2/b1*delta_t*mue2pix^2;
% this leads to an incremental decay value of
decay_val=exp(-delta_t/tau);



%How many files do we have, this sets the index for the fow loop below
max_s=length(flow_files);


% Create the folder for the results!
%save_path=uigetdir(path,'Where should I store the results??');
mkdir(path,'\internal_forces');
mkdir(path,'\internal_forces_im');


%Now we need to estimate how long we should take the relaxing stresses into
%account. We will drop them, once the values are smaller than 1%
ignore_small=ceil(log(0.01)/log(decay_val))+1;

%now we create the array that will hold the stress history. This needs to
%be as big as ux and uy, and the 3rd dimension (which represents time in this case) needs to be "ignore_small"
%long. To get the right size, we are just reading once the ux and so on.
load([path,'\',flow_files(3).name]);
ux=x_flow;
uy=y_flow;
[s_y,s_x]=size(ux);
Fx(1:s_y,1:s_x,1:ignore_small)=0;
Fy(1:s_y,1:s_x,1:ignore_small)=0;
i_s=ignore_small;

for i=1:max_s
    tic
    %Here we load the displacement fields
    load([path,'\',flow_files(i).name]);
    ux=x_flow;
    uy=y_flow;
    %Since in the elastic force calculation, we do everything in pixel
    %units, we need to transform the ux, and uy (which were saved in µm/min
    %back into pixel between 2 frames.
    %Thus we have to multiply with delta_t/60*1/mue2pix
    ux=ux*delta_t/60*1/mue2pix;
    uy=uy*delta_t/60*1/mue2pix;
    %now calculate the i'th stress

    ux(find(ux==0))=NaN;
    uy(find(uy==0))=NaN;
    [fx_int,fy_int]=lin_elast_force(ux,uy,E,g,mue2pix);
    
    Fx(:,:,i_s)=fx_int;
    Fy(:,:,i_s)=fy_int;
    
    
    
    %memorize the actual sum stress, and save in the "viscoleastic_stress
    %folder
    sig_x=sum(Fx,3);
    sig_y=sum(Fy,3);
    
    save([path,'\internal_forces\internal_force_array',num2str(i+1000-1),'.mat'],'sig_x','sig_y');

    
    %Now generate the force image and save in internal_forces_im
    %First plot the radial force part
    [sig_th,sig_r]=cart2pol(sig_x,sig_y);
    h=figure;
    imshow(sig_r)
    colormap('jet')
    set(gca,'CLim',[0,Max_F])
    colorbar
    hold on

    %now get the flow field arrows
    plot_data=generate_plot_normalized(sig_x,sig_y,flow_field_arrow_distance);
    quiver(plot_data(:,1),plot_data(:,2),plot_data(:,3),plot_data(:,4),0.3,'y');
    print('-dpng','-r300',[path,'\internal_forces_im\','display_vis_',num2str(i+1000-1),'.png']);
    %saveas(h,[path,'\internal_forces\',num2str(index),'.fig']);
    hold off
    close
    
    
 %   save_in_image(sig_x,sig_y,[path,'\internal_forces'],i);
    
    %shift the stress according to deformation
    
    Fx_sh=shift_according_to_deformation_3d(ux,uy,Fx);
    Fy_sh=shift_according_to_deformation_3d(ux,uy,Fy);
    
    % let the gel relax according to decay_val    
    Fx_sh=decay_val*Fx_sh;
    Fy_sh=decay_val*Fy_sh;
    
    Fx=Fx_sh;
    Fy=Fy_sh;
    
    %and finally shift the Fx, and Fy arrays in the 3rd DImension, to get
    %rig of the oldest value and make space at the i_s position
    
    Fx(:,:,i_s+1)=0;
    Fx(:,:,i_s+1)=0;
    
    Fx(:,:,1:i_s)=Fx(:,:,2:i_s+1);
    Fx(:,:,1:i_s)=Fx(:,:,2:i_s+1); 
    

    
    toc
    pause(2)
end

    
    