f = fopen('C:\MATLAB\LIB\browse\CD.txt');
CD = fread(f);
fclose(f);

cd(char(CD'));
