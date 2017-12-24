CD = cd;
f = fopen('C:\MATLAB\LIB\browse\CD.txt','w');
fwrite(f,CD);
fclose(f);

cd('c:\DATA\DataTest');