function displayComp(Pts,discComp,useImage,sizeIm)
  
  FALSE = (0 == 1);
  TRUE = ~FALSE;

  colourSize = 33; 
  fcMap = hsv(colourSize);  
  leafComp = discComp; 
  figure; clf; 
  if (useImage) 

     mx = 0.5;

     im = Pts(:,1);
     rgb0 = repmat(im, [1 3])/255;

     rgb = mx * rgb0;

     for k=1:size(leafComp,2) 
       idx = leafComp(:,k); 
       if any(idx) 
         c = fcMap(1+mod(floor(rand*95357),colourSize),:);
         rgb(idx,:) = rgb(idx,:) + ...
                      (1-mx) * rgb0(idx, :) .* repmat(c, [sum(idx) 1]);
       end
     end

     rgb = reshape(rgb, [sizeIm 3]);
     image(rgb);
     axis off; axis equal;

  else
      
    qPts = Pts(:, [1,2]); 
    plot(Pts(:,1), Pts(:,2),'.r'); 
    hold on; 
    for k=1:size(leafComp,2) 
      idx = leafComp(:,k); 
      if any(idx) 
        c = 1+mod(floor(rand*95357),colourSize); 
      
        ht= text(qPts(idx,1), qPts(idx>0,2),sprintf('%d',k), ... 
	 'Color', fcMap(c,:), 'FontName', 'Courier', 'FontWeight', 'bold', ... 
	       'FontSize', 20); 
      end
    end 
  end 
  hold off; 
  
 
