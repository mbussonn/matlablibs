function [F_x,F_y]=lin_elast_force(ux,uy,E,g,mue2pix)
%[F_x,F_y]=lin_elast_force(ux,uy)
%uses lineae elasticity theory to estimate the forces that deform an
%elastic gel
%first get the strain tensor field from the deformation fields
%then calculate the stress tensor field following linear elasticity theory
%and Landau-Lifshitz. Finally we use the equilibrium condition to identify
%the sum of internal and external forces from the spatial derivative of the
%stress tensor
%We do all the calculation in pixel units, and as the last step, we
%transform back into micrometer. Thus ux, and uy should be given in pixel
%units!

%mode to double
ux=double(ux);
uy=double(uy);

%first check if E and g are given, if not use defaults E=1000, g=0.4
if (nargin < 3)
    delta_t=1;
    pixtomue=1;
    E=1000;
    g=0.4;
end


!

% use a 3D array, to display the strain and stress tensor, Elements are
% represented as follows:
% 11 12
% 21 22
%means:  strain(:,:,1)=11;strain(:,:,2)=12;strain(:,:,3)=21;strain(:,:,4)=22;
[ux_x,ux_y]=gradient(ux);
[uy_x,uy_y]=gradient(uy);
strain(:,:,1)=ux_x;
strain(:,:,2)=0.5*(ux_y+uy_x);
strain(:,:,3)=strain(:,:,2);
strain(:,:,4)=uy_y;


%Now use linear Elasticity theory from Landai Lifshitz to get the stress
%tensors
stress(:,:,1)=E/((1+g)*(1-2*g))*((1-g)*strain(:,:,1)+g*strain(:,:,4));
stress(:,:,4)=E/((1+g)*(1-2*g))*((1-g)*strain(:,:,4)+g*strain(:,:,1));
stress(:,:,3)=E/(1+g)*strain(:,:,3);
stress(:,:,3)=stress(:,:,2);

%now use the equilibriom condition to establish the force-fields
%d stress_ik/d x_k=-F_i
[grad_str_xx_x,grad_str_xx_y]=gradient(stress(:,:,1));
[grad_str_xy_x,grad_str_xy_y]=gradient(stress(:,:,2));
[grad_str_yy_x,grad_str_yy_y]=gradient(stress(:,:,4));

F_x=(-grad_str_xx_x-grad_str_xy_x)./mue2pix.^2;
F_y=(-grad_str_xy_y-grad_str_yy_y)./mue2pix.^2;

