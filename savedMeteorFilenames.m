%{
This script will give the saved meteor filenames.  Be sure to update it if
you make more saved meteors!
%}

disp('Creating the meteor filenames from savedMeteorFilenames...')

opticalFilenames = {'0330105633ultra.mat',
    '0330105754ultra.mat',
    '0330110015ultra.mat',
    '0330111636ultra.mat',
    '0330112439ultra.mat',
    '0330114821ultra.mat',
    '0330121306ultra.mat',
    '0330122329ultra.mat',
    '0330123903ultra.mat', %start the 31st data next line
    '0331131516ultra.mat',
    '033112598ultra.mat',
    '033112589ultra.mat',
    '0331124154ultra.mat',
    '0331124035ultra.mat',
    '0331122239ultra.mat',
    '0331121857ultra.mat',
    '0331113559ultra.mat',
    '0331112933ultra.mat',
    '0331112056ultra.mat',
    '0331104642ultra.mat',
    '033110388ultra.mat',
    '0331101723ultra.mat',
    '033110244ultra.mat',
    '033195531ultra.mat',
    '033194147ultra.mat',
    '033193226ultra.mat',
    '033192643ultra.mat',
    '033192558ultra.mat',
    '03319722ultra.mat',
    '033182858ultra.mat'};

radarFilenames = {'0331105634.mat',
    '0331105754.mat',
    '033111015.mat',
    '0331111637.mat',
    '0331112440.mat',
    '0331114822.mat',
    '033112137.mat',
    '0331122329.mat',
    '033112394.mat', %start of 31st data next line
    '0331131516.mat',
    '033112598.mat',
    '033112589.mat',
    '0331124154.mat',
    '0331124035.mat',
    '0331122239.mat',
    '0331121857.mat',
    '0331113559.mat',
    '0331112933.mat',
    '0331112056.mat',
    '0331104642.mat',
    '033110388.mat',
    '0331101723.mat',
    '033110244.mat',
    '033195531.mat',
    '033194147.mat',
    '033193226.mat',
    '033192643.mat',
    '033192558.mat',
    '03319722.mat',
    '033182858.mat'};

% useful information for the meteors
totalFramesAllMeteors = [13, 26, 14, 34, 55, 25, 41, 19, 17, 25];
startFramesAllMeteors = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20,1,];
goodHoughIndeces = {[2 3 4 5];
    [2 3 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ];
    [ ]};