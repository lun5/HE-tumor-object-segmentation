function [] = examineRGB(rgbIm)
%
% Plot an rgb image, do a PCA of the colour channels, and plot the
% the results.
% Input: an rgb image, stored in an ny x nx x 3 double array.  The
% values in rgb must be in [0,1].

  if size(size(rgbIm),2) ~= 3
    fprintf(1, 'Expecting a nx x ny x 3 double array.');
    return;
  end
  if max(max(max(rgbIm))) > 1.0 | min(min(min(rgbIm))) < 0.0
    fprintf(1, 'Expecting rgb values in [0,1].');
    return;
  end
  
  sizeIm = size(rgbIm(:,:,1));

  % Display the colour image (figure 1).
  figure(1); clf;
  image(rgbIm); 
  axis('image'); axis('off');
  % display  at double size: 2 screen pixels/image pixel
  if max(size(rgbIm)) < 250
    resizeImageFig(figure(1), size(rgbIm(:,:,1)), 2);  
  else
    resizeImageFig(figure(1), size(rgbIm(:,:,1)), 1); 
  end 
  hold off;
  drawnow;

  % Display the individual r, g and b channels.
  % The gray level range in each channel is reduced by 1/sqrt(3) so that
  % they can be compared directly with the PCA images below.
  figure(2);  clf;
  grayMap = gray(256);
  colormap(grayMap);
  channel = {'Red'; 'Green'; 'Blue'};
  for k=1:3
    h = subplot(2,2,k);
    image(rgbIm(:,:,k)*255/sqrt(3));  % Scale down to match PCA image.
    axis('image'); axis('off');
    title(channel(k));
  end
  hold off;
  drawnow;
  % Note the 3 images look quite similar, indicating that there is a
  % significant amount of correlation between the r,g, and b channels.
  % The differences between the channels is apparent only in the
  % strongly coloured regions.

  % Investigate the correlations by doing a PCA in the rgb space.
  rgb = reshape(rgbIm, prod(sizeIm), 3);
  mRgb = mean(rgb);
  rgb = rgb - repmat(mRgb, size(rgb, 1), 1);
  C = (rgb' * rgb)/prod(sizeIm);
  [U, S, V] = svd(C, 0);
  if ones(3,1)' * U(:,1) < 0
    U(:,1) = -U(:,1);
    V(:,1) = -V(:,1);
  end 
  var = diag(S);
  sigma = sqrt(var);
  fprintf(1, 'sigma = %6.4f %6.4f %6.4f\n', sigma);

  % Transform the rgb components to the principal directions
  rgb = rgb * U;
  compMax = max(max(abs(rgb)));
  compMax = max(compMax, 1);

  %% The set of rgb values looks like a 3D elliptical pancake or cigar...
  %% Show the scatter plot of the rgb data projected
  %% onto the second and third components, versus the first component.
  figure(3); clf;
  skip = 38;
  for cDim = 2:3
    subplot(1,2,cDim-1);
    plot(rgb(1:skip:prod(sizeIm), 1), rgb(1:skip:prod(sizeIm), cDim), '.r'); 
    axis([-compMax compMax -compMax compMax]);
    title(sprintf('rgb components 1 and %d', cDim));
    xlabel('Coeff. of 1st comp.');
    ylabel(sprintf('Coeff for comp. # %d', cDim));
    hold on;
    theta = [-1:0.01:1]*pi;
    x = sigma(1)*cos(theta);
    y = sigma(cDim) * sin(theta);
    plot(x,y,'g', 2*x, 2*y, 'b');
    axis square;
    hold off;
  end
  drawnow;
  %% The curves show 1 and 2 standard deviation ellipses for
  %% the Gaussian model.  The Gaussian model is sometimes a poor approximation
  %% of the data, especially for the more constrained scenes of just
  %% a few objects.  For these images the scatter plots look clumpy.
  %% The clumps correspond to common colours in the image.

  % Display the principal component channels.
  % The gray level range in each channel is reduced by 1/sqrt(3) to
  % allow for the maximum possible value.  This is the same scale factor
  % as for the rgb images in figure 2, so they can be compared directly.
  rgb = rgb + repmat(mRgb, size(rgb,1), 1);
  rgb = reshape(rgb, [sizeIm, 3]);
  figure(4); clf;
  grayMap = gray(256);
  colormap(grayMap);
  channel = {'Red'; 'Green'; 'Blue'};
  for k=1:3
    h = subplot(2,2,k);
    image(rgb(:,:,k)*255/sqrt(3)); % Scale down to avoid saturation.
    axis('image'); axis('off');
    str = sprintf('Comp#%d: %4.2fR', k, U(1,k));
    for j=2:3
      if (U(j,k) < 0.0)
        str = [str sprintf('%5.2f%s', U(j,k), channel{j}(1))];
      else
        str = [str sprintf('+%4.2f%s', U(j,k), channel{j}(1))];
      end
    end
    title(str);
  end
  drawnow;
  % The PCA has arranged that the correlation between these
  % components is zero.
  % Note that the range of the first component is larger than the
  % range of the individual r, g, and b channels, while the range of
  % the second and third components are smaller.  
  % Also, notice that the first component is brightness, some weighted
  % combination of R, G, and B, with positive weights. The second component
  % typically has opposite signs for the R and B weights (eg. R or Y vs B),
  % while the third component has the same sign weights for R and B and
  % the opposite sign for G (eg. R or P (purple) vs G).

  clear rgb U V C;
