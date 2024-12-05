% Written by Jeremiah V John

function saveToVideo(V, fname)
%% Saves given 3D array to a greyscale MP4 file of chosen name
%       V (uint8): 3D array of dimensions (pixel, pixel, frame)
%       fname (string [or char list whatever Matlab does]): chosen name

writer = VideoWriter(fname, 'MPEG-4');
open(writer);
for i=1:size(V,3)
    writeVideo(writer,mat2gray(V(:,:,i)));
end
close(writer);
end



