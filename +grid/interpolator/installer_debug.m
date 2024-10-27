function installer_debug(installfolder)

    cd(installfolder);
    % make object file
    mex  interpUtil.cpp -c -g
    % link
    mex -g ppmval.cpp interpUtil.obj 
    % compiling ppuval
    mex -g ppuval.cpp interpUtil.obj
    
end