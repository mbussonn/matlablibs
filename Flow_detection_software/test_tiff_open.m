fname = 'E:\Science\data\vasily\CompositeTruncated.tif';
info = imfinfo(fname);
num_images = numel(info);
for k = 1:num_images
    A = imread(fname, k, 'Info', info);
    % ... Do something with image A ...
end