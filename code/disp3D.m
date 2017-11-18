clear
%% read stacked frame files
Frames = 100;
StackNum = 38;

ratFrm = 0.5;

if exist('stack.mat')
    load stack;
else
    for k=1:Frames*ratFrm
        for l=1:StackNum

            DigitDiff=floor(log10(Frames))-floor(log10(k));   
            %fprintf('Frame:%d Stack:%d \n',k,l)
            
            if k == 1
             DigitDiff=floor(log10(Frames));
            end

            if DigitDiff == 0
            Stack(:,:,l,k)=uint16(imread(['stack_' int2str(k) '.tif'],l));
            end

            if DigitDiff == 1
            Stack(:,:,l,k)=uint16(imread(['stack_0' int2str(k) '.tif'],l));
            end

            if DigitDiff == 2
            Stack(:,:,l,k)=uint16(imread(['stack_00' int2str(k) '.tif'],l));
            end

            if DigitDiff == 3
            Stack(:,:,l,k)=uint16(imread(['stack_000' int2str(k) '.tif'],l));
            end
        end
        h = waitbar(k/Frames/ratFrm);
    end
end
    
save('stack.mat','Stack','-v7.3')
close(h);
%Stack=double(Stack);
Stack=Stack-min(Stack(:));
sliceomatic(uint8(Stack(:,:,:,1)));