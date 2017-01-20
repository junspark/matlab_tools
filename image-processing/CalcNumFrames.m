function num_frame = CalcNumFrames(bytes, buffer_size, frame_size)
num_frame   = (bytes - buffer_size)/frame_size;
return