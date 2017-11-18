

if ~exist('inputInfo.mat')
    fname4D = input('Enter the name of the ORIGINAL file you wanna analyze (with the extension)   ', 's');
        
    %READ Number of Frames
    Frames = input('Enter the number of Frames');

    %READ Window Size
    WindowSize = input('Enter the Window Size (pixels)   ');
    BigWindowSize=WindowSize+4;

    %READ Pixel Size
    PixelSize= input('Enter the Pixel Size (nanometers)   ');

    %READ Number of Stacks
    StackNum = input('Enter the number of Z-Planes   ');

    %READ Distance between Planes
    PlaneDist = input('Enter the distance between Planes (nanometers)   ');
    
    save('inputInfo.mat','WindowSize', 'Frames', 'PixelSize', 'StackNum', 'PlaneDist', 'fname4D')
else 
    load inputInfo;
end


zt2stacks;
TraCKer_3D_w_ZcolorPlot_deep_strain;
shearDisp;