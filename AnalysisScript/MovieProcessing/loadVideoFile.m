function V = loadVideoFile(fname)
%% Loads a video file, converts to grayscale and returns the 3D vector (HxWxT)
    v = VideoReader(fname);
    V = zeros(v.Height, v.Width, v.NumFrames, 'uint8');
    
    for i = 1:v.NumFrames
        V(:,:,i) = rgb2gray(read(v,i));
    end
end