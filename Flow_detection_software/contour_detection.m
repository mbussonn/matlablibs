function contour_detection
%This functin tries to detect the countour of an aggregate.
%First we blur the image, then we try to correct for inhomegneious
%illumination, then we try to dtermine a smart threshold, and finally we
%create a binary mask that will be used in the main program

size_g=55;%size of the gaussian kernel
sig=1;%sigma of the gaussian kernel
d=10;%edge rechtalnges to compensate inhomogenous illumination
t=1;%should I invert?
ones1=5;%first imclose, removes smal particle
ones2=50;%second imclose, removes particles at the end of procedure and clese the contour
min_size=50;%remove particles after threshold

%here I load the image
im=imread('E:\Science\data\sarah_geraldo\FasG3 sph phase\FasG3 sph phase0000.png');

%I want that center is bright, so inverse if necessary
if t==1
    im=double(im);
    im=uint16((im-2^16)*-1);
end

%Now I correct for assymtric illumination. I take the corners and
%intermolate a map for normalisation

[x,y]=size(im);
I(1,1)=mean(reshape(im(1:d,1:d),1,d^2));
I(2,2)=mean(reshape(im(x-d+1:x,y-d+1:y),1,d^2));
I(2,1)=mean(reshape(im(1:d,y-d+1:y),1,d^2));
I(1,2)=mean(reshape(im(x-d+1:x,1:d),1,d^2));

%now I interpolare the normalisation matrix
[X,Y] = meshgrid(1:x-1:x,1:y-1:y);
Z = I;
[XI,YI] = meshgrid(1:x,1:y);
ZI = interp2(X,Y,Z,XI,YI);

%now I normalize the image
im_n=im-(uint16(ZI')-min(min(ZI)));

%# Create the gaussian filter with hsize = [5 5] and sigma = 2
G = fspecial('gaussian',[size_g size_g],sig);
%# Filter it
im_blur = imfilter(im_n,G,'same');
% now threshold
graythresh(im_blur)
im_th=im2bw(im_blur, graythresh(im_blur));

%now I clean up
bw2 = imfill(im_th,'holes');
bw3 = imopen(bw2, ones(ones1,ones1));
bw4 = bwareaopen(bw3, min_size);
bw4b=(imclose(bw4,ones(ones2,ones2)));
bw5=(imfill(bw4b,'holes'));

%finally remove all particles that don't take up half of the total white
%pixels numbers
bw6= bwareaopen(bw5, round(sum(sum(bw5))/2));

bw6_perim = bwperim(bw6);
overlay1 = imoverlay(im, bw6_perim, [1 0 0]);
imshow(overlay1)

%# Display

imshow(overlay1)

