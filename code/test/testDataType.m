clear all;
N=23; K = 1e5; Nz=floor(N/2); % zero element param.
gDbl= double(rand(N,N,N,N)*K);
ZEROgDbl = gDbl;
ZEROgDbl(1:Nz,:,:,:) = 0;
save('gDbl.mat','gDbl');
save('ZEROgDbl.mat','ZEROgDbl');
whos; clear gDbl ZEROgDbl
gSngl= single(rand(N,N,N,N)*K);
ZEROgSngl = gSngl;
ZEROgSngl(1:Nz,:,:,:) = 0;
save('gSngl.mat','gSngl');
save('ZEROgSngl.mat','ZEROgSngl');
whos; clear gSngl ZEROgSngl
gLgc = logical(rand(N,N,N,N)*K);
ZEROgLgc = gLgc;
ZEROgLgc(1:Nz,:,:,:) = 0;
save('gLgc.mat','gLgc');
save('ZEROgLgc.mat','ZEROgLgc');
whos; clear gLgc ZEROgLgc
gUINT8 = uint8(rand(N,N,N,N)*K);
ZEROgUINT8 = gUINT8;
ZEROgUINT8(1:Nz,:,:,:) = 0;
save('gUINT8.mat','gUINT8');
save('ZEROgUINT8.mat','ZEROgUINT8');
whos; clear gUINT8 ZEROgUINT8
gUINT16 = uint16(rand(N,N,N,N)*K);
ZEROgUINT16 = gUINT16;
ZEROgUINT16(1:Nz,:,:,:) = 0;
save('gUINT16.mat','gUINT16');
save('ZEROgUINT16.mat','ZEROgUINT16');
whos; clear gUINT16 ZEROgUINT16

%% 
for p = 1:1e3
    a = 2^p+1;
    aDouble = double(a);
    aSingle = single(a);
    if rem(aDouble,100) ~= rem(aSingle,100)
        p
        break
    end
end
 % p =24
 
 %uint8 => < 256
 %uint16 => < 65536
 