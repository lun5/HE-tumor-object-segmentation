vers = version;
mx = mexext;
if ispc
  mx = 'dll';  % Use dll extension for both Matlab 6 and 7
end

eval(['mex corrDn.c wrap.c convolve.c edges.c -output corrDn.' mx]);
eval(['mex upConv.c wrap.c convolve.c edges.c -output upConv.' mx]);
eval(['mex pointOp.c -output pointOp.' mx]);
eval(['mex histo.c -output histo.' mx]);
eval(['mex range2.c -output range2.' mx]);
