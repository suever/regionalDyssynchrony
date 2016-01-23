function compile_mex()
    % compile_mex - Compiles all regional dyssynchrony mex-files

    % Copyright (c) 2013, Jonathan Suever
    % Last Modified: 02-15-2013
    % Modified By: Suever (suever@gmail.com)

    basedir = fileparts(mfilename('fullpath'));

    mex('-largeArrayDims', ...
        fullfile(basedir, 'mex', 'xcorrPeriodic.c'), ...
        '-lmwblas', '-O', '-outdir', basedir);

    mex('-largeArrayDims', ...
        fullfile(basedir, 'mex', 'qtClustering.c'), ...
        '-O', '-outdir', basedir);
end
