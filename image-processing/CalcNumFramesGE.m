function num_frames = CalcNumFramesGE(pfname)

num_frames  = dir(pfname);
num_frames  = (num_frames.bytes-8192)/2048/2048/2;