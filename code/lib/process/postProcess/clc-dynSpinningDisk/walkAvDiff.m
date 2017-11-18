% calculates the intensity diff. data from the acquired movies
clear all
fnacqDIR = rdir('clc\acq_*');
fnacq = fnacqDIR.name;
fnacq2DIR = rdir('dyn\acq_*');
fnacq2 = fnacq2DIR.name;

ix = strfind(fnacq2,'acq');

fnacqout = [fnacq(1:ix+2) 'Diff' fnacq(ix+3:end)];
fnacq2out = [fnacq2(1:ix+2) 'Diff' fnacq2(ix+3:end)];
delete(fnacqout);
delete(fnacq2out);

imginfo = imfinfo(fnacq);
frm1=numel(imginfo);
imginfo = imfinfo(fnacq2);
frm2=numel(imginfo);

for i = 1:frm1 % read acquisition images
    A_ = imread(fnacq,i); % clathrin channel
    A(:,:,i) = A_;
end

for i = 1:frm2 % read acquisition images
    D_ = imread(fnacq2,i); % dynamin channel
    D(:,:,i) = D_;
end



img0 = A; 
WAsz = 3;
%% walking average code
nEl = size(img0,1)*size(img0,2);
img=reshape(img0,nEl,size(img0,3));
WA = double(ones(1,WAsz));
WA = WA/sum(WA);
imgWA = conv2(double(img),WA,'valid');
imgWA=reshape(imgWA,size(img0,1),size(img0,2),size(img0,3)-WAsz+1);
%%
Awa = imgWA;
clear imgWA;


img0 = D; 
WAsz = 2;
%% walking average code
nEl = size(img0,1)*size(img0,2);
img=reshape(img0,nEl,size(img0,3));
WA = double(ones(1,WAsz));
WA = WA/sum(WA);
imgWA = conv2(double(img),WA,'valid');
imgWA=reshape(imgWA,size(img0,1),size(img0,2),size(img0,3)-WAsz+1);
%%
Dwa = imgWA;


%% difference images
Adiff = Awa(:,:,1:end-4) - Awa(:,:,5:end); % CLC  (+ + + _ - - - )
Ddiff = Dwa(:,:,1:end-2) - Dwa(:,:,3:end); % DYN  ( + + - -)

%% dbg
if 0
            t0 = t-4; if t0<1, t0=1; end;
            t0_2 = t-1; if t0_2<1, t0_2=1; end;
            
            t1 = t+1; if t1>frm2, t1=frm2; end;
            t2 = t+4; if t2>frm2, t2=frm2; end;

            A0 = mean(A(:,:,t0:t0_2),3);
            A0_2 = A0;
            A2 = mean(A(:,:,t1:t2),3);
            A02temp = A0-A2;
            A02temp(A02temp<0)=0;
            figure(1345);imagesc(A02temp); axis image
end
%%

%% save diff images
for i = 1:size(Adiff,3)
    imwrite(uint16(Adiff(:,:,i)),fnacqout,'WriteMode','append');
end

for i = 1:size(Ddiff,3)
    imwrite(uint16(Ddiff(:,:,i)),fnacq2out,'WriteMode','append');
end

