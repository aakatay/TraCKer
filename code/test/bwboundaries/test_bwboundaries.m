  %% Example 0
  
  a = zeros(10);
  a(2:8,2:8)=0.1;
  %a(2:5,2:5)=1;
  a(4:6,4:6)=1;
  a(5,5)=1.1;
  a(9,9)=0.1;
  aBW = imregionalmax(a,4); 
  figure(1);imagesc(aBW);
  B=bwboundaries(aBW);
  aBW = imregionalmax(a);
  numel(B)
  B{1}
  figure(2);
  imagesc(a); hold on
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end  
    hold off
    
    %% Example 1
    
%         Read in and threshold the rice.png image. Display the labeled
%         objects using the jet colormap, on a gray background, with region
%         boundaries outlined in white.
 
        figure(1)
       I = imread('rice.png');
       BW = im2bw(I, graythresh(I));
       [B,L] = bwboundaries(BW,'noholes');
       imshow(label2rgb(L, @jet, [.5 .5 .5]))
       hold on
       for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
       end
       hold off
 
    %% Example 2
    
%     Read in and display binary image blobs.png. Overlay the region 
%     boundaries on the image. Display text showing the region number
%     (based on the label matrix), next to every boundary. Additionally,
%     display the adjacency matrix using SPY.
 
%     HINT: After the image is displayed, use the zoom tool in order to read
%           individual labels.
    figure(2)
       BW = imread('blobs.png');
       [B,L,N,A] = bwboundaries(BW);
       imshow(BW); hold on;
       colors=['b' 'g' 'r' 'c' 'm' 'y'];
       for k=1:length(B),
         boundary = B{k};
         cidx = mod(k,length(colors))+1;
         plot(boundary(:,2), boundary(:,1), colors(cidx),'LineWidth',2);
         %randomize text position for better visibility
         rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
         col = boundary(rndRow,2); row = boundary(rndRow,1);
         h = text(col+1, row-1, num2str(L(row,col)));
         set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
       end
       figure(3); spy(A);
 
    %% Example 3
    
%     Display object boundaries in red and hole boundaries in green.
 figure(4)
       BW = imread('blobs.png');
       [B,L,N] = bwboundaries(BW);
       imshow(BW); hold on;
       for k=1:length(B),
         boundary = B{k};
         if(k > N)
           plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
         else
           plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
         end
       end
       
    %% Example 4
    
%     Display parent boundaries in red (any empty row of adjacency
%     matrix belongs to a parent) and their holes in green.
 figure(5)
       BW = imread('blobs.png');
       [B,L,N,A] = bwboundaries(BW);
       imshow(BW); hold on;
       for k=1:length(B),
         if(~sum(A(k,:)))
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
           for l=find(A(:,k))'
             boundary = B{l};
             plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
           end
         end
       end