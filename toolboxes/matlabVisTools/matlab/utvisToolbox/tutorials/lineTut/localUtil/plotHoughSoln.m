function [xLine, yLine] = plotHoughSoln(im, corners, H, aExtent, thetaDeg, x0, minVotes)

tH = H;
tH(tH(:,:) < minVotes) = 0;
tH(tH(:,:) >= minVotes) = 1;
a = zeros(sum(sum(tH)), 2);

[y x] = size(H);
boxX = (aExtent(2)-aExtent(1))/x;
boxY = (aExtent(4)-aExtent(3))/y;
k = 1;
for i=1:y
  for j=1:x
    if (H(i,j) > minVotes)
      a(k,1) = aExtent(1) + (j-0.5)*boxX;
      a(k,2) = aExtent(3) + (i-0.5)*boxY;      
      k = k + 1;  
    end
  end
end

thetaRad = thetaDeg * pi/180;
nrml = [cos(thetaRad) sin(thetaRad)];
tang = [-nrml(2) nrml(1)];

[x y] = size(a);
n0 = repmat(nrml, x, 1);
t0 = repmat(tang, x, 1);

p = n0 + repmat(a(:,2), 1, 2) .* t0;
n = p ./ norm(p);
c = a(:,1) ./ norm(p);

yLine = repmat([corners(1,2) corners(2,2)], x, 1);
xLine = [ ((-n(:,2) .* (yLine(:,1) - x0(2)) - c(:,1)) ./ n(:,1)  + x0(1)) ...
          ((-n(:,2) .* (yLine(:,2) - x0(2)) - c(:,1)) ./ n(:,1)  + x0(1))];
hold on;
for i=1:1:x
 line(xLine(i,:), yLine(i,:));
end
hold off;
