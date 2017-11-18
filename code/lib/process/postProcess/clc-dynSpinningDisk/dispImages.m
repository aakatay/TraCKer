% disp images for profile selection (in selectStruct)


        tc=fc(i); % c-chan peak time
        td=fd(i); % d-chan peak time
        t = td;
        x = xxc(i);
        y = yyc(i);
        if ~isDisp & ismember(tc,skipFrames), isContinue=1; return; end;
        if exist('xytLim')
            if x<xytLim(1) | x>xytLim(3) | y<xytLim(2) | y>xytLim(4), isContinue=1; return; end;
        end
        
        ix1d = find(fd == td);
        ix2d = find(fd == td+1);
        ix1c = find(fc == td);
        ix2c = find(fc == td+1);
        xd1All = xxd0(ix1d);
        yd1All = yyd0(ix1d);
        xc1All = xxc0(ix1c);
        yc1All = yyc0(ix1c);
        xd2All = xxd0(ix2d);
        yd2All = yyd0(ix2d);
        xc2All = xxc0(ix2c);
        yc2All = yyc0(ix2c);
        
        if round(x)<2 | round(y)<2 | round(x)>size(Cimg,2)-2 | round(y)>size(Cimg,1)-2 % skip at edges
            profC0(i,:)=nan;
            profC2(i,:)=nan;
            profClate(i,:)=nan;
            profD0(i,:)=nan;
            profD2(i,:)=nan;
            profDlate(i,:)=nan;
            isContinue=1;return;
        end

        if tDiffSingleFrame
            t0 = t-8; if t0<1, t0=1; end;
            t0_2 = t-1; if t0_2<1, t0_2=1; end;
            t1 = t+1; if t1>frm2, t1=frm2; end;
            t2 = t+2; if t2>frm2, t2=frm2; end;
            t3 = t+3; if t3>frm2, t3=frm2; end;
            tLate = t3+20; if tLate>frm2, tLate=frm2; end;

            % clathrin intensity averages at diff. times around dyn peak
            A0 = mean(A(:,:,t0:t),3);
            A0_2 = mean(A(:,:,t0_2:t),3);
            A2 = mean(A(:,:,t1:t1),3);
            A3 = mean(A(:,:,t2:t2),3);
            A4 = mean(A(:,:,t3:t3),3);
            Alate = mean(A(:,:,t3:tLate),3);
            D0 = mean(D(:,:,t0:t),3);
            D0_2 = mean(D(:,:,t0_2:t),3);
            D2 = mean(D(:,:,t1:t1),3);
            D3 = mean(D(:,:,t2:t2),3);
            D4 = mean(D(:,:,t3:t3),3);
            Dlate = mean(D(:,:,t3:tLate),3);  
        else % walking average data
            %%
            dbg = 0;
            if dbg, t=4; end;
            t0 = t-3; if t0<1, t0=1; end;
            t0_2 = t-1; if t0_2<1, t0_2=1; end;
            
            t1 = t+1; if t1>frm2, t1=frm2; end;
            t2 = t+3; if t2>frm2, t2=frm2; end;
            
            t3 = t+3; if t3>frm2, t3=frm2; end;
            tLate = t3+20; if tLate>frm2, tLate=frm2; end;

            % clathrin intensity averages at diff. times around dyn peak
            A0 = mean(A(:,:,t0:t0_2),3);
            A0_2 = A0;
            A2 = mean(A(:,:,t1:t2),3);
            A3 = nan(size(A2));
            A4 = A3;
            Alate = mean(A(:,:,t3:tLate),3);
            %%
            
            D0 = mean(D(:,:,t-1:t),3);
            D0_2 = D0;
            D2 = mean(D(:,:,t+1:t+2),3);
            D3 = nan(size(D2));
            D4 = D3;
            Dlate = mean(D(:,:,t3:tLate),3);     
            if dbg
                A02temp = A0-A2;
                A02temp(A02temp<0)=0;
                figure(1345);imagesc(A02temp); axis image
                %set(gca,'Xlim',xl)
                %set(gca,'Ylim',yl)
                colormap('gray')
                
                D02temp = D0-D2;
                D02temp(D02temp<0)=0;
                figure(1346);imagesc(D02temp); axis image
                %set(gca,'Xlim',xl)
                %set(gca,'Ylim',yl)
                colormap('gray')
            end     
        end

        %if y<40, continue; end; % false recs
        xc0 = xxc0(i); % double precision
        yc0 = yyc0(i);
        xd0 = xxd0(i); % double precision
        yd0 = yyd0(i);
        x0 = xc0-fs; if x0<1, x0=1; end
        x1 = xc0+fs; if x1>size(Cimg,2), x1=size(Cimg,2); end;
        xl = [x0 x1];
        y0 = yc0-fs; if y0<1, y0=1; end
        y1 = yc0+fs; if y1>size(Cimg,1), y1=size(Cimg,1); end;
        yl = [y0 y1];
        
        x0 = round(xc0-fs2); if x0<1, x0=1; end
        x1 = round(xc0+fs2); if x1>size(Cimg,2), x1=size(Cimg,2); end;
        xl2 = [x0 x1];
        y0 = round(yc0-fs2); if y0<1, y0=1; end
        y1 = round(yc0+fs2); if y1>size(Cimg,1), y1=size(Cimg,1); end;
        yl2 = [y0 y1];

if isLoadProfData & size(exy,1)>=i 
    if ~isReAdjustProfLines, isContinue=1;return; end
            cx=exy(i,1);
            cy=exy(i,2);
            ex=exy(i,3);
            ey=exy(i,4);
            if isnan(exy(i,3))
                cx = xd0 - 1; cy = yd0;
                ex = (xd0-cx)*2+cx; % endpoints
                ey = (yd0-cy)*2+cy;
            end
            %if isnan(ex), isContinue=1;return; end;
        else
            cx = xd0 - 1; cy = yd0;
            ex = (xd0-cx)*2+cx; % endpoints
            ey = (yd0-cy)*2+cy;
        end
        e1x=cx; e2x=ex; e1y=cy; e2y=ey;

  

        %tvec = t0:t2+5;

        %% images
        if ~isSelectStructs
            figure(97) % clathrin images before and after dyn peak
            amx = max([A0(:); A2(:); A3(:)]);
            dmx = max([D0(:); D2(:); D3(:)]);
            colormap(parula(256))
            sc = 256/amx;
            sc2 = 256/dmx;
            subplot(4,3,1)
            image(A0*sc); axis image; title('A0(t0:t)'); 
            subplot(4,3,2)
            image(A2*sc); axis image; title('A2(t+1)'); 
            subplot(4,3,3)
            image(A3*sc); axis image; title('A3(t+2)'); 
            subplot(4,3,4)
            imagesc(D0*sc2); axis image; title('D0(t0:t)'); 
            subplot(4,3,5)
            imagesc(D2*sc2); axis image; title('D2(t+1)'); 
            subplot(4,3,6)
            imagesc(D3*sc2); axis image; title('D3(t+2)'); 

            subplot(4,3,7)
            imagesc(Alate); axis image; title('Alate(t3:t3+20)'); 
            subplot(4,3,8)
            imagesc((A0-A2)); axis image; title('A0-A2'); colorbar
            subplot(4,3,9)
            imagesc((A2-A3)); axis image; title('Adiff2(t+2)-(t+1)'); colorbar
            subplot(4,3,10)
            imagesc(Dlate); axis image; title('Dlate(t3:t3+20)'); 
            subplot(4,3,11)
            D02 = D0-D2; D02(D02<0)=0;
            imagesc((D02)); axis image; title('D0-D2');  colorbar
            subplot(4,3,12)
            D23 = D2-D3; D23(D23<0)=0;
            imagesc((D23)); axis image; title('Ddiff2(t+2)-(t+1)'); colorbar
        
        else % isSelectStructs
        %
            Asel = double(A(yl2(1):yl2(2),xl2(1):xl2(2),t-3:t+3));
            Dsel = double(D(yl2(1):yl2(2),xl2(1):xl2(2),t-3:t+3));
            Asel = uint16(Asel/max(Asel(:))*256);
            Dsel = uint16(Dsel/max(Dsel(:))*256);
            
            figure(977)
            CM = parula(256);
            colormap(CM)
            subplot(7,2,1)
            image(Asel(:,:,1)); axis image;
            subplot(7,2,3)
            image(Asel(:,:,2)); axis image;
            subplot(7,2,5)
            image(Asel(:,:,3)); axis image;
            subplot(7,2,7)
            image(Asel(:,:,4)); axis image;
            subplot(7,2,9)
            image(Asel(:,:,5)); axis image;
            subplot(7,2,11)
            image(Asel(:,:,6)); axis image;
            subplot(7,2,13)
            image(Asel(:,:,7)); axis image;
            
            subplot(7,2,2)
            image(Dsel(:,:,1)); axis image;
            subplot(7,2,4)
            image(Dsel(:,:,2)); axis image;
            subplot(7,2,6)
            image(Dsel(:,:,3)); axis image;
            subplot(7,2,8)
            image(Dsel(:,:,4)); axis image; 
            subplot(7,2,10)
            image(Dsel(:,:,5)); axis image; 
            subplot(7,2,12)
            image(Dsel(:,:,6)); axis image; 
            subplot(7,2,14)
            image(Dsel(:,:,7)); axis image; 
            btn = uicontrol('Style', 'pushbutton', 'String', 'print',...
            'Position', [20 20 50 20],...
            'Callback', 'printSelStructs');  
        
            Aint = Asel(4:6,4:6,:);
            aint = mean(reshape(Aint,9,7),1);
            Dint = Dsel(4:6,4:6,:);
            dint = mean(reshape(Dint,9,7),1);
            figure(988)
            plot([aint;dint]'); legend('clc','dyn')
            
        end

if isDisp, return; end;
        %%
        %clf(100)
        delete(100)
        figure(100)
        imagesc(A0); axis image;

        % xy line profile
        [px1,py1,p0_1] = improfile(A0_2,[e1x xd0],[e1y yd0]);
        [px2,py2,p0_2] = improfile(A0_2,[xd0 e2x],[yd0 e2y]);