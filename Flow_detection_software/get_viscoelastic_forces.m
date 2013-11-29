function Pol=get_viscoelastic_forces(path)
%using the displacement fields, we calculate the sresses that are required
%to produce these deformtions. This is done by assuming a stress relaxation
%function that can be either modeled by a exponential decay, or be measured
%with an AFM. The viscoelastic responce is the nothing but the decay
%function multiplied with the Elasticity coefficient in the of the pure
%elastic case, multiplied this the strain rate, integrated over time!

% Okay, we do this the folowing way:
% First read the first deformations, and calculate the according stresses (in Array S),
% then shift them, and multiply with exp(-delta_t/tau).
% After that, do the next images, and put in it in S(:,:,2), the sum along
% the 3rd dimension is the resulting stress!
% Afterthat repeat the shift and the relaxation. Do this over and over
% untill done. 

%First we need to get the rigth path, if not given at call
%if nargin < 1, 
 %   path=uigetdir('D:\','Give me the basepath with the Info file of the data to analyze.');
%end
path='E:\backup_leipzig\E\Retro_flow_for_thesis\GFP_Act_NG108_050204 fresh_control_4 [++ Stat 600sec]\retrograde_flow_data'

%path='D:\Data_retr\GFP-Actin NG108 050205 Big Stationary [done,++]\GFP-Actin NG108 050205 Big Stationary [done,++]\standard\retrograde_flow_data';%If you change this, change the mue2pix ratio too


%Here we define the coefficients we got out of the AFM-measurements for the
%growth cone, using the dashpot - spring||dashpot model. These values are
%measured with YunBe. They correspond to: E=112, tau=2.4, eta=4.8
a1=905.7536;
a2=17.6628;
b1=8.5519


%Then we try to extract some important parameters:
%First of all, try to read the Info_GUI file, to get pixle to µm ratio and the time between frames:
info_gui_list=dir([path,'\Info_GUI_v*.txt'])
if length(info_gui_list)==1
    fid=fopen([path,'\',info_gui_list(1).name],'rt')
    y = 0;
    while feof(fid) == 0
       tline = fgetl(fid);
      matches = findstr(tline, 'Voxel size in the images in µm:');
      matches_time=findstr(tline, 'Time between frames in sec:');
      num = length(matches);
       if num > 0
          [token,rem] = strtok(tline,':');
          mue2pix=str2double(rem(2:length(rem)));
       end
      num_time = length(matches_time);
       if num_time > 0
           [token,rem] = strtok(tline,':');
          %time between image frames
          delta_t=str2double(rem(2:length(rem)));
       end           
    end
    fclose(fid); 
else
    error('Sorry, but I could not find the Info_GUI file. Please create it, and store the µm to pix ration in it')
end


%Set some important parameters:



%According to the dashpot- spring||dashpot model, the prefactor before the
%exponential decay reads:
E=(a1*b1-a2)/b1^2;
%NOw we convert it to pixel units:
E=E*mue2pix^2;
%Poisson Ratio
g=0.47;


%relaxation time
tau=b1;
%since the measured deformation allway tkes place during the delta_t time
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
x_flow_im_name_in=dir([path,'\dir_flow_analysis\x*.png']);
max_s=length(x_flow_im_name_in);


% Create the folder for the results!
%save_path=uigetdir(path,'Where should I store the results??');
mkdir(path,'\internal_forces');


%Now we need to estimate how long we should take the relaxing stresses into
%account. We will drop them, once the values are smaller than 5%

ignore_small=ceil(log(0.01)/log(decay_val))+1;

%now we create the array that will hold the stress history. This needs to
%be as big as ux nad uy, and the 3rd Dimension needs to be "ignore_small"
%long. To get the right size, we are just reading once the ux and so on.
[ux,uy,I]=arrange_and_display_folder(path,1);
[s_y,s_x]=size(ux);
Fx(1:s_y,1:s_x,1:ignore_small)=0;
Fy(1:s_y,1:s_x,1:ignore_small)=0;
i_s=ignore_small;

for i=1:max_s
    tic
    %Here we load the actual image, and the displacement fields
    [ux,uy,I]=arrange_and_display_folder(path,i);
    %Since in the elastic force calculation, we do everything in pixel
    %units, we need to transform the ux, and uy (which were saved in µm/min
    %back into pixel between 2 frames.
    %Thus we have to multiply with delta_t/60*1/mue2pix
    ux=ux*delta_t/60*1/mue2pix;
    uy=uy*delta_t/60*1/mue2pix;
    %now calculate the i'th stress
    [fx_int,fy_int]=lin_elast_force(ux,uy,E,g,mue2pix);
    
    Fx(:,:,i_s)=fx_int;
    Fy(:,:,i_s)=fy_int;
    
    
    
    %memorize the actual sum stress, and save in the "viscoleastic_stress
    %folder
    sig_x=sum(Fx,3);
    sig_y=sum(Fy,3);
    
    save_in_image(sig_x,sig_y,[path,'\internal_forces'],i);
    
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

    
    