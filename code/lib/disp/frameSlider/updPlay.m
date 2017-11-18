
playSpeed = get(hplaySpeedSlider,'value'); % fps
playState = get(hPlay,'String');
isPlaying = 1;
isPlayWin = 0;
framePlayTime = 0;
df1 = 0; dfnum = 0;
if strcmp('playState','play')
    set(hPlay,'String','stop');
    framePlayTime = 1/playSpeed;
    dfnum = 1;
elseif strcmp('playState','play window')
    set(hPlay,'String','stop');
    framePlayTime = 1/playSpeed;
    df1 = 1;
    isPlayWin = 1;
elseif strcmp('playState','stop')
    set(hPlay,'String','play');
end