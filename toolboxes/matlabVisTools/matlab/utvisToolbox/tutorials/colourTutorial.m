%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File: colourTutorial.m
%%%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This Matlab tutorial examines the spectral composition of diffusely
% reflected light given natural light sources and materials.  It introduces
% the CIE xy-colour space and the XYZ colour matching functions
% used to characterize perceptually equivalent chromatic sources.
% The CIE (i.e. the International Commission on Illumination (en Francais),
% http://www.cie.co.at/ciecb/), provides the dominant international
% standards on illumination.
%
% The relationship of the CIE colour matching functions to the
% spectral sensitivites of the colour sensors in the human retina
% is described.  We also show that the set of spectral reflectances
% for 'natural' materials under standard daylights
% can be well approximated within a three dimensional space.  This
% motivates why our retinas have just three types of colour sensors
% with different spectral properties.   
% 
% Finally, the presence of correlations in natural reflectance distributions,
% and in the sensor responses, motivates looking for similar correlations
% between the R, G, and B channels of digital images.  We illustrate
% the presence of such correlations, and briefly discuss their properties.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prerequisite: The singular value decomposition.  Eg. see the
%%%               svdTutorial from the Standford iseToolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Dependencies: 
%       iseToolbox/pyrTools/
%          histo.m
%       iseToolbox/color/
%          chromaticity.m   data/  daylightSpectrum.m
%       iseToolbox/color/data/
%          D65.mat	  daycie.mat             hit489.mat
%          sanyoLCD.mat	  smithPokornyCones.mat  xyz.mat
%
%       utvisToolbox/colour/
%          colourTutorial.m    xy2xyImage.m
%          parseMunsellName.m  examineRGB.m
%       utvisToolbox/colour/data/
%          munsellNew.mat  spectrum.tif  diagxy.tif
%       utvisToolbox/file
%          resizeImageFig.m
%       images/tif    %% You can replace these by your own tif images.
%          hedvig3_050.tif  allanA_080.tif
%          hats.tiff        lego-microserf.tif
%          im7h.tif
%
%%%  Author: ADJ, Fall 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;

%% Set the root for the vision toolbox.
%% This is already done if you ran the matlabVisRoot startup file.
global matlabVisRoot
if size(matlabVisRoot,1) == 0
  matlabVisRoot = '~jepson/pub/matlab';
end

%% Set the root of the different toolbox directories.
iseRoot = [matlabVisRoot '/iseToolbox'];
utvisRoot = [matlabVisRoot '/utvisToolbox'];
imageRoot = [matlabVisRoot '/images'];

%% Set the path to include these directories
path(path, [iseRoot '/pyrTools/MEX']);
path(path, [iseRoot '/pyrTools']);
path(path, [iseRoot '/color']);
path(path, [iseRoot '/color/data']);
path(path, [utvisRoot '/colour']);
path(path, [utvisRoot '/colour/data']);
path(path, [imageRoot '/tif']);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab data sets for natural light sourses and the Munsell
% chip reflectances.  These are explored further below.
%%%%%%%%%%%%%%%%%%%%%%%%%%
load([iseRoot '/color/data/xyz.mat']);
load([utvisRoot '/colour/data/munsellNew.mat']);
load([iseRoot '/color/data/D65.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VISIBLE SPECTRUM
%%%%%%%%%%%%%%%%%%%%%%%%%%
% The visible spectrum corresponds to wavelengths from 400 to 700 nm
% roughly.  The data in this tutorial will be represented over
% a slightly larger range [370, 730], sampled every 10 nm. 
mnW = min(wavelength);
mxW = max(wavelength);
display([mnW mxW]);

%% A picture of the visible spectrum 
%% (thanks to the Colour Group: http://www.city.ac.uk/colourgroup/)
figure(1); clf;
[im map] = imread('data/spectrum.tif');
image(im); colormap(map);
% Label axes and put on ticks
set(gca, 'YTickLabel', '');
set(gca, 'XTick', [13:60:193]);
set(gca, 'XTickLabel', [400:100:700]);
xlabel('(Apprx) Wavelength (nm)');
title('Visible Spectrum');
drawnow;
clear im map;
% The colours in this spectrum image are NOT calibrated, 
% since the display is not calibrated.  But this image gives
% a rough indication of the association between various
% wavelengths and their colours.  A more accurate association
% is provided in the CIE plots considered further below.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STANDARD DAYLIGHTS
%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CIE (http://www.cie.co.at/ciecb/) has constructed a model for
% standard daylights.  One typical choice is
% daylight D65.  Here the 65 denotes the associated
% colour temperature of 6500 degrees Kelvin.
plot(wavelength, D65);
title('Spectral Power Distribution of Standard Daylight D65');
axis([mnW mxW 0 200]);
drawnow;

% Other daylights are specified by their   
% colour temperature, tDeg in [4, 25]K degrees.
% lower temperatures -> 'redder' illuminants near sunrise and sunset.
% higher temperatures -> 'bluer' illuminants on clear sunny days.
% All the daylights have been normalized to 100 at wavelength 560nm.
tDeg=[4:2:10];
figure(1); clf;
plot(wavelength, D65, 'k'); hold on;
axis([350 750 0 200]);
for k=1:size(tDeg,2)  
 daylight = daylightSpectrum(tDeg(k) * 1000);
 plot(wavelength, daylight, 'r');
 axis([mnW mxW 0 200]);
 title(sprintf('SPD tDeg = 4-%dK (r), 6.5K (k)', tDeg(k)));
 pause(1);  
end

%% We use the standard daylight D65 in the following computations.
%% D65 is the same as the standard daylight with colour temperature 6500
tDeg = 6.5;
daylight = daylightSpectrum(tDeg * 1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REFLECTANCES: MUNSELL CHIPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Reference: Munsell Book of Color - Matte Finish Collection 
%%  (Munsell Color, Baltimore, Md., 1976). 
%% The Munsell set consists of a large collection of callibrated
%% colour chips (as in paint 'chips'), whose colours were systematically 
%% chosen to span a wide range of possibilities.  The number of chips
%% in the data set is
size(munsell, 2)

%% Each column of the matrix munsell provides the samples
%% of the spectral reflectance distribution of one Munsell chip
%% for the range of wavelengths given in the vector wavelength.
%% Here wavelength ranges from:
fprintf(1, 'wavelength range: %dnm to %dnm, in steps of %dnm\n',mnW,mxW, ... 
        wavelength(2) - wavelength(1));

%% Lets take a look at some of these reflectance distributions.
% Plot every 50th reflectance
figure(1); clf;
skip = 50;
skipInd = [1:skip:size(munsell,2)];
plot(wavelength, munsell(:,skipInd));
axis([mnW mxW 0 1.1]);
xlabel('Wavelength (nm)');
ylabel('Reflectivity');
title('Various Munsell Reflectances');
drawnow;
%% These reflectance distributions appear to have been extrapolated
%% by a constant value between wavelengths 370 and 380.  This approximation
%% will not significantly effect the calculations below.

%% The Munsell chips are named for their perceived colour:
for k=skipInd
   fprintf(1,'k: %d, name: %s\n', k, munsellNames(k,:)); 
end

%% These Munsell names describe the hue, brightness and saturation
%% for the particular chip.
help parseMunsellName;

%% For example, we can select the indices for munsell reflectances in
%% the series '5G' (medium green, with various lightnesses and saturation)
selectInd = [];
for k = 1:size(munsell,2)
  [hueMod hueName lght sat] = parseMunsellName(munsellNames(k,:));
  if strcmp(hueName, 'G') & hueMod == 5
    selectInd = [ selectInd , k];
  end       
end

% echo the selected names  (note they are all in the series 5G
% with various lightnesses an saturations
munsellNames(selectInd,:)

%
% Then plot the selected reflectance distributions.  
% Here you need to type any key in the Matlab window to see the 
% next reflectance in the series.  
%
figure(1); clf;
cnt = 0;
for k = selectInd
  plot(wavelength, munsell(:,k));
  title(['Munsell Chip: ' munsellNames(k, :)]);
  axis([370 730 0 1.1]);
  xlabel('Wavelength (nm)');
  ylabel('Reflectivity');
  if cnt < 5
    % fprintf(1, 'Type any key to continue...\n');
    cnt = cnt+1;
  end
  pause(1);     
  % pause(0.2);
end
fprintf(1, '...done.\n');
% Note the plot titles indicate that these selected indices all have
% hue '5G', as requested above.
% Also, notice the full munsell names are of the form: 
%    5G lght/sat
% where lght corresponds to the overall lightness of the
% reflectance and sat corresponds to the saturation of the
% colour (i.e. how much different from white it is).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STIMULI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Illuminating these chips with a natural daylight gives a
%% spectral distribution function of 
%%     s(lambda) = r(lambda) * I(lambda)
%% Here we take the term L dot N to be constant (i.e. a flat
%% surface) and ignore this constant term.

stimuli = munsell .* repmat(daylight, [1 size(munsell, 2)]);

%% Plot every 50th reflectance
figure(1); clf;
skip = 50;
skipInd = [1:skip:size(munsell,2)];
plot(wavelength, stimuli(:,skipInd));
xlabel('Wavelength (nm)');
ylabel('Power');
title('Spectral Power for Munsell Reflectances in Daylight');
drawnow;
%% Try changing the daylights, and replotting the reflectances.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SENSOR SPECTRAL SENSITIVITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([iseRoot '/color/data/smithPokornyCones.mat']);
%% spCones is a 361 x 3 matrix providing the spectral sensitivities
%% of the 3 types of cones (colour sensors) for an average human.
%% The peak values of the spectral sensistivities have been scaled
%% to be one.  Variations in the transmittance of the lens
%% alters the spectral sensitivities across individuals, and for
%% the same individual at different ages.
figure(1); clf;
plot(wavelength,spCones(:,3),'b', wavelength, spCones(:,2),'g',...
     wavelength, spCones(:,1),'r');
title('Cones Spectral Sensitivities');
xlabel('Wavelength(nm)');
ylabel('Sensitivity');
drawnow;
%% We see that there are three flavours of cones, the short (S), 
%% medium (M) and long wavelength (L).

% Compute the cone responses to each of the stimuli.  This
% is just the integral of the product of the cone sensitivity times
% the spectral power over wavelengths:
resp = spCones' * stimuli;

% Show a scatter plot of responses
figure(1); clf;
subplot(1, 2, 1);
plot(resp(1,:), resp(2,:), 'ob'); axis square;
axis([0 10000 0 10000]);
title('Linear Cone Responses: L vs M');
xlabel('L-cone response');
ylabel('M-cone response');
subplot(1, 2, 2);
plot(resp(1,:), resp(3,:), 'ob'); axis square;
axis([0 10000 0 10000]);
title('Linear Cone Responses: L vs S');
xlabel('L-cone response');
ylabel('S-cone response');
drawnow;

% Note that the L- and M-cone responses are much more correlated than
% the M- and S-cone responses.  The discrete steps that are clear
% in the L vs M plot are due to the Munsell reflectances having roughly a
% discrete set of mean reflectances, i.e. the predominantly integer
% lightnesses between 2 and 10.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2+1D COLOUR SPACE and METAMERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It would be useful to characterize a given stimulus in terms of the
% power absorbed by the L, M, and S cones.  Stimuli which give
% identical powers for the three cones should then be perceptually
% identical. A pair of such stimuli are called metamers.  We will
% discuss this more further below.
 
% For now, note that the spectral powers absorbed by the
% cones are just three numbers, and that various values of these three
% numbers would then characterize the effect of all of the stimuli we
% have considered so far.  Moreover, if we were to normalize
% for the brightness somehow, say by dividing by the
% sum of the absorbed power, then we would get just two
% numbers to characterize the 'colour' of a stimulus.  Any pair of
% stimuli which provided the same numbers, that is, the
% same relative power in each of the sensors, would be metamers.

% In other words, we can precisely represent each of the stimuli...
stimuli = munsell .* repmat(daylight, [1 size(munsell, 2)]);
figure(1); clf; 
subplot(2,2,1);
skip = 50;
skipInd = [1:skip:size(munsell,2)];
plot(wavelength, stimuli(:,skipInd)); 
axis([mnW mxW 0 100]); axis square;
xlabel('Wavelength (nm)');
ylabel('Power');
title('Munsell * D65');
drawnow;
%... by the absorbed power for the three cones (3D), or by
% the brightness scaled cone absorbances (2D):
resp = spCones' * stimuli;
totPow = sum(resp,1);
for k=skipInd
 fprintf(1, '%s: %6.1f %6.1f %6.1f: %4.2f %4.2f %4.2f\n', ...
   munsellNames(k,:), resp(:,k)', resp(:,k)'/totPow(k));
end

% For example, after scaling for brightness, these same Munsell stimuli
% can be characterized by the relative power absorbed by the three cone
% types.  This replaces each spectral distribution by just one point (!)
% in each of the following scatter plots:
figure(1); subplot(2,2,3);
plot(resp(1,skipInd)./totPow(skipInd), resp(2,skipInd)./totPow(skipInd), 'ob');
axis square; axis([0 1 0 1]);
title('Fraction of Total Power');
xlabel('L-cone');
ylabel('M-cone');
subplot(2,2,4);
plot(resp(1,skipInd)./totPow(skipInd), resp(3,skipInd)./totPow(skipInd),'ob');
axis square; axis([0 1 0 1]);
title('Fraction of Total Power');
xlabel('L-cone');
ylabel('S-cone');
hold off;
drawnow;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CIE XYZ-COLOUR COORDINATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We stated above that if we adjust two light sources with different
% spectral distribution functions such that the spectral power
% absorbed by the three cone types are the same for
% both light sources, then they are said to be
% 'metamers'.  We cannot distinguish metamers.  The consequence
% of this is that we can describe the chromaticity of lights in terms
% of only 2 parameters (or 3 if brightness is included).

% The CIE has carefully calibrated 2 and 3D colour spaces
% in order to specify perceptually identical chromatic stimuli.  
% There are several different CIE spaces that are all closely
% related.  Here we will use the CIE XYZ colour matching functions,
% and the related CIE xy-space.  Later below we will show that these colour
% matching functions are closely related to the spectral sensitivities
% of human cones, as expected.

% To convert any spectral distribution to CIE XYZ coords,
% we need to use the XYZ-colour matching functions.  These
% are simply three functions of wavelength. 
% We can plot them:
figure(1); clf;
plot(wavelength,XYZ(:,1),'b', wavelength, XYZ(:,2),'g',...
     wavelength, XYZ(:,3),'r');
title('XYZ-color matching functions');
xlabel('Wavelength(nm)');
drawnow;
 
% To compute the XYZ coordinates of a spectral distribution
% we need to integrate the products of these colour matching functions
% with the given distribution over the visible spectrum.  In Matlab
% we can approximate these three integrals by simply multiplying the
% distribution by the matrix XYZ'

% Computing the XYZ coords of all of the munsell reflectances
% illuminated by the standard daylight:
XYZcoords = XYZ' * stimuli;

% Each column of this product gives the X,Y,Z coordinates
% for the spectral distribution, eg:
XYZcoords(:,1)

% The importance of the XYZ coordinates is that, for the average
% human observer, stimuli with the same XYZ coordinates are
% perceptually identical.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  CIE xy-COLOUR COORDINATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CIE x,y coordinates are formed by dividing the X and Y coords
% by the sum of the X,Y, and Z coords.  This scales out the overall
% brightness.
xyMunsCoords = XYZcoords(1:2,:) ./ repmat(sum(XYZcoords, 1), 2, 1);

% Eg. the first column has XYZ coords:
XYZcoords(:,1)'
% with the sum
sum(XYZcoords(:,1))
% and  x,y coordinates:
xyMunsCoords(:,1)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CIE xy-COLOUR SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here is a colour image of the xy-colour space is:
%% (thanks to the Colour Group: http://www.city.ac.uk/colourgroup/)
figure(1); clf;
[im map] = imread('data/diagxy.tif');
image(im); colormap(map); axis off;
title('CIE x,y Colour Space');
drawnow;
% The colours in this image are NOT calibrated, 
% since the display is not calibrated.  But this image gives
% a rough indication of the association between various
% xy-colour coordinates and their perceived colours.

% This 2D x,y space does not indicate the overall brightness
% of a light source, but rather only the hue and saturation.
% The overall brightness is included in the related XYZ coordinates.  

% We can superimpose the xy-colour coordinates we computed for
% the Munsell set (illuminated by D65) on the colour image for
% the xy-space.  
% To do this we need to convert CIE xy-coords to the corresponding pixel in
% the diagxy.tif image.  We have implemented this coordinate
% conversion in the function xy2xyImage: 
xyMunsImage = xy2xyImage(xyMunsCoords);

% And then a scatter plot can be done using these pixel coordinates. 
% Assuming selectInd has been computed as above, giving the
% indices of the '5G' hue series in the Munsell set...
figure(1); hold on;
scatter(xyMunsImage(1,selectInd), xyMunsImage(2,selectInd), 5, 'w');
hold off;
drawnow;
%% The white circles indicate the CIE xy coords of the selected
%% Munsell chips.

%% Try selecting different subsets of the Munsell set
%% and replotting.  For example, to select only the
%% series with hueMod = 5 (no matter what the hueName is). 
selectInd = [];
for k = 1:size(munsell,2)
  [hueMod hueName lght sat] = parseMunsellName(munsellNames(k,:));
  if hueMod == 5
    selectInd = [ selectInd , k];
  end       
end
figure(1); clf;
[im map] = imread('data/diagxy.tif');
image(im); colormap(map); axis off;
title('CIE x,y Colour Space');
hold on;
scatter(xyMunsImage(1,selectInd), xyMunsImage(2,selectInd), 5, 'w');
hold off;
drawnow;
%% You should see ten radial arms for each of R, YR, Y, GY, G, BG, B,
%% PB, P, RP.
%% Try variations on this, such as requiring the saturation value
%% sat >= 8.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHANGING THE ILLUMINANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the whole Munsell set, illuminated by the standard daylight D65
%% (tDeg = 6.5) we get
figure(1); clf;
[im map] = imread('data/diagxy.tif');
image(im); colormap(map); axis off;
title('CIE x,y Colour Space');
hold on;
scatter(xyMunsImage(1,:), xyMunsImage(2,:), 5, 'w');
hold off;
drawnow;

%% Changing the daylight to sunset/sunrise:
tDeg = 4.0;
redDaylight = daylightSpectrum(tDeg * 1000);
% Repeating the computation of the XYZ coords of all the stimuli...
tmpXYZ = XYZ' * (munsell .* repmat(redDaylight, [1 size(munsell, 2)]));
% ... converting XYZ -> xy...
tmpxy = chromaticity(tmpXYZ);
% ... and then xy -> pixel coords for the colour xy-space image...
tmpxyImage = xy2xyImage(tmpxy);
% ... and finally plotting the points:
figure(1);
hold on;
scatter(tmpxyImage(1,:), tmpxyImage(2,:), 5, 'r');
hold off;
drawnow;
% Notice how the redder daylight shifted the xy coords of all of the stimuli
% upwards and to the right.

%% Change the daylight to a bright blue sky:
tDeg = 20.0;
blueDaylight = daylightSpectrum(tDeg * 1000);
tmpXYZ = XYZ' * (munsell .* repmat(blueDaylight, [1 size(munsell, 2)]));
tmpxy = chromaticity(tmpXYZ);
tmpxyImage = xy2xyImage(tmpxy);
hold on;
scatter(tmpxyImage(1,:), tmpxyImage(2,:), 5, 'b');
hold off;
drawnow;
% Notice how the blue-ish daylight shifted the xy coords of all
% of the stimuli downwards and to the left.
clear tmp*;

%% Colour constancy refers to our ability to infer the
%% reflectance properties of objects (eg. identify which Munsell chip)
%% and neglect (or account for) the variability of the stimulus
%% due to various light sources having different spectral intensity
%% distributions.  That is, colour constancy somehow allows us to correct
%% for this shift in the xy coords due to changing illuminates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WHERE DO THESE CIE XYZ-COLOUR MATCHING FUNCTIONS COME FROM?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The essential mathematical property of the three XYZ-colour matching
% functions is that they span the same linear space as the three
% the 'effective' spectral sensitivities of the human cones.  Here,
% 'effective' refers to the spectral sensitivities INCLUDING 
% the absorbances due to the lens and macular pigment in the eye.

% The spectral sensitivities of individual cone cells are determined
% by genetics, and are the same across humans with normal colour vision.  
% However, different people have different absorption spectrums for
% light absorbed by the eye itself before it reaches the cones (eg.
% by the lens and macular pigment).  Unfortunately for us theoreticians
% and engineers, this eye absorption spectrum changes with each individual,
% changes with age, and even changes across the retina within a given
% individual.  However, the variability is small enough for CIE to get a
% useful approximation with one standard observer.

% So what are the effective (after eye absorbtion) spectral sensitivities
% for this CIE standard observer?   

% We can show that the Smith-Pokorny spectral sensitivities,
% which include the effects of eye absorption, are approximately 
% linear combinations of XYZ
size(XYZ)
[Uxyz S Vxyz] = svd(XYZ, 0);
Sxyz = diag(S);
XYZpseudo = (Vxyz * S^-1) * Uxyz';  % see the svdTutorial for pseudo-inv

figure(1);clf;
plot(wavelength,spCones(:,3),'b', wavelength, spCones(:,2),'g',...
     wavelength, spCones(:,1),'r');
title('Smith-Pokorny cone senistivites');
xlabel('Wavelength (nm)');
ylabel('Sensitivity');
axis([mnW  mxW 0 1.1]);
% fprintf(1, 'Type any key to continue...\n');
pause(1);

%% Use the pseudo-inverse of XYZ to compute the least squares
%% fit of the spCones...
spXYZcoeffs = XYZpseudo * spCones;
xyzSensors = XYZ * spXYZcoeffs;
hold on;
plot(wavelength,xyzSensors(:,3),'c', wavelength, xyzSensors(:,2),'y',...
     wavelength, xyzSensors(:,1),'m');
title('Smith-Pokorny cone senistivites, fit to XYZ');
hold off;
drawnow;
% Note that the Long and Medium wavelength cone sensitivities are
% approximated very well by linear combinations of the XYZ basis
% vectors.  The approximation is less accurate for the
% Short wavelength cone.  This is due to different assumptions
% about the lens and macular pigment absorptions in the `average
% eye'.  In the following we will use the CIE XYZ average observer
% as the standard observer, and use these xyzSensors as the spectral
% sensitivities for this average observer.

%% This approximation is negligible in terms of the
%% relative sensitivities to natural stimuli.
resp = spCones' * stimuli;
resp2 = xyzSensors' * stimuli;
% Show a scatter plot of responses
chLabel='LMS';
figure(1); clf;
for k=1:3
  subplot(2, 2, k);
  plot(resp(k,:), resp(k,:), '.b'); axis square;
  axis([0 10000 0 10000]);
  xlabel(sprintf('SP Channel %c', chLabel(k)));
  ylabel(sprintf('xyz Channel %c', chLabel(k)));
end
drawnow
% The approximation of using xyzSensors in place of Smith-Pokorny
% sensors makes the math relating XYZ to cone power absorptions
% exact.  For example, stimuli with the same XYZ
% colour coordinates have EXACTLY the same effects on
% the xyzSensors.  This relation would only be approximate
% if we used the Smith-Pokorny sensors.  In effect, the xyzSensors
% are effective spectral sensitivities that could have corresponded
% exactly to the 'average observer' used to generate the XYZ-colour
% matching functions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARTS OF THE CIE xy-COLOUR SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the CIE xy-space white/gray is indicated at (x,y) = (0.33, 0.33). 
% For example, here are some grey stimuli (uniform spectral distribution)
grayStim = ones(size(wavelength)) * [0.2: 0.2: 1.0];
figure(1); clf;
plot(wavelength, grayStim);
title('Ideal White and Gray Reflectances');
xlabel('Wavelength(nm)');
ylabel('Reflectivity');
axis([370 730 0 1.1]);
hold off;
drawnow;
% And here are the XYZ coords (one column for each brightness)...
XYZgrayStim = XYZ' * grayStim
% ... along with the corresponding xy coords.
xyStim = chromaticity( XYZ' * grayStim )

% Lines moving radially outward from the white/gray point 
% in xy-space correspond to a particular hue. 
% The saturation increases with distance to the white point. 
% Therefore, as we saw above, the '5G' series corresponds to points
% (roughly) on a radial line starting at the white point.

% The curving portion of the boundary of the CIE xy space corresponds 
% to monochromatic stimuli (i.e. the visible spectrum, or
% laser light with a particular wavelength).  The straight
% line across the bottom-right corresponds to mixtures of
% two monochromatic stimuli at the extremes of the visible
% spectrum, (roughly) 370nm and  730nm.
% To demonstrate these facts, construct some mono-chromatic distributions:
selectWave = [400:50:700]';
% Choose the array indices k so that wavelength(k) in array selectWave.
selectK = selectWave-wavelength(1)+1;  
% Check that we did the selection correctly.
wavelength(selectK)'
%%... or, alternatively ...
all(wavelength(selectK) == selectWave)
% Build the monochromatic stimuli
stimuli = eye(size(wavelength,1));  %% identity matrix
% The k-th column of stimuli has a one in the kth row, zero
% elsewhere.  This corresponds to a laser light with wavelength
% wavelength(k).

% Plot the selected monochromatic stimuli
figure(1); clf;
for k = selectK'
  figure(1);
  plot(wavelength, stimuli(:,k));
  title(sprintf('Monochromatic Stimulus: %d',wavelength(k)));
  axis([370 730 0 1.1]);
  xlabel('Wavelength (nm)');
  ylabel('Reflectivity');
  pause(0.5);
end

% Compute xy coords for ALL of these monochromatic stimuli
xyStimuli = chromaticity( XYZ' * stimuli);

%% Plot the xy coords for these laser stimuli, and place
%% dots at the wavelengths 400:50:700, as specified by selectWave.
figure(2); clf; %% Separate the windows for figures 1 and 2 on your screen
[im map] = imread('data/diagxy.tif');
image(im); colormap(map); axis off;
title('CIE x,y Colour Space');
hold on;
xyImage = xy2xyImage(xyStimuli);
plot(xyImage(1,:), xyImage(2,:), 'w');
scatter(xyImage(1,selectK), xyImage(2, selectK), 5, 'w');
hold off;
drawnow;
%% Here the white curve corresponds to the xy coords of the monochromatic
%% sources.  The dots on this curve correspond to wavelengths 
wavelength(selectK)'
%% ordered with 400nm (violet/blue) at the bottom left,
%% moving clockwise with increasing wavelength to 700nm (red)
%% at the bottom right corner.

% The straight line across the bottom of the CIE xy space
% corresponds to mixtures between two extreme monochromatic sources.
redLaser = zeros(size(wavelength));
redLaser(size(wavelength,1)) = 1.0;
redLaser = redLaser / sum(XYZ' * redLaser); 
%% The division above scales X+Y+Z = 1.  This is used to
%% match the intensity of the light to other lights used in the mixture.
violetLaser = zeros(size(wavelength)); 
violetLaser(1) = 1.0;
violetLaser = violetLaser / sum(XYZ' * violetLaser); %% Scales X+Y+Z =1
mix = 0.0:0.2:1.0;
stimuli = zeros(size(wavelength,1), size(mix,2));
for k = 1:size(mix,2)
  stimuli(:,k) = mix(k) * redLaser + (1.0-mix(k)) * violetLaser;
end  
%% Peak at the mixture stimuli
figure(1); clf;
for k = 1:size(mix,2)
  figure(1);
  plot(wavelength, stimuli(:,k));
  title(sprintf('Mixtures of Monochromatic Stimuli: mix = %0.1f',mix(k)));
  axis([350 750 0 500]);
  xlabel('Wavelength (nm)');
  ylabel('Reflectivity');
  pause(0.5);     
  % pause(1);
end
% Add these points to the x,y-space plot.  They form
% the sloped straight line on the bottom of the plot.
figure(2); hold on;
xyStimuli = chromaticity(XYZ' * stimuli);
xyImage = xy2xyImage(xyStimuli);
plot(xyImage(1,:), xyImage(2,:), 'w');
scatter(xyImage(1,:), xyImage(2, :), 5, 'w');
hold off;
drawnow;

% Similarly, we can make a range of synthetic reddish lights 
% having different saturations as follows.  First make a white
% light having a flat spectral distribution.
whiteLight = ones(size(wavelength));
%% Scale the light so the XYZ coords sum to one.  This
%% is the same scaling used for the redLaser above.
whiteLight = whiteLight / sum(XYZ' * whiteLight); 
mix = 0.0:0.2:1.0;
stimuli = zeros(size(wavelength,1), size(mix,2));
%% Build the mixture stimuli
for k = 1:size(mix,2)
  stimuli(:,k) = mix(k) * redLaser + (1.0-mix(k)) * whiteLight;
end  
%% Take a peak at these stimuli
figure(1); clf;
for k = 1:size(mix,2)
  figure(1);
  plot(wavelength, stimuli(:,k));
  title(sprintf('Red Monochromatic + White Stimuli: mix = %0.1f',mix(k)));
  axis([350 750 0 0.02]);
  xlabel('Wavelength (nm)');
  ylabel('Reflectivity');
  pause(0.5);
  % pause(1);
end
% Add these points to the x,y-space plot:
figure(2); hold on;
xyStimuli = chromaticity(XYZ' * stimuli);
xyImage = xy2xyImage(xyStimuli);
plot(xyImage(1,:), xyImage(2,:), 'w');
scatter(xyImage(1,:), xyImage(2, :), 5, 'w');
hold off;
drawnow;

% Similarly, we can make a range of synthetic greenish lights 
% having various different saturations as follows:
greenLaser = zeros(size(wavelength));
greenLaser(510-wavelength(1)+1) = 1.0; 
greenLaser = greenLaser / sum(XYZ' * greenLaser);
mix = 0.0:0.2:1.0;
stimuli = zeros(size(wavelength,1), size(mix,2));
for k = 1:size(mix,2)
  stimuli(:,k) = mix(k) * greenLaser + (1.0-mix(k)) * whiteLight;
end  
%% Take a peak at these stimuli
figure(1); clf;
for k = 1:size(mix,2)
  figure(1);
  plot(wavelength, stimuli(:,k));
  title(sprintf('Green Monochromatic + White Stimuli: mix = %0.1f',mix(k)));
  axis([350 750 0 0.02]);
  xlabel('Wavelength (nm)');
  ylabel('Reflectivity');
  pause(0.5);
  % pause(1);
end
%% Note the reflectivity of the uniform (white) component
%% in this mixture must be much lower at each wavelength
%% than the green laser is at it's peak in order for the
%% two components to have the same X+Y+Z sum.  So in these
%% plots we have truncated the y-axis to show the variation
%% of the white component only.
%
% Add these points to the x,y-space plot:
figure(2); hold on;
xyStimuli = chromaticity(XYZ' * stimuli);
xyImage = xy2xyImage(xyStimuli);
plot(xyImage(1,:), xyImage(2,:), 'w');
scatter(xyImage(1,:), xyImage(2, :), 5, 'w');
hold off;
drawnow;

%% Super-imposing the '5G' munsell line formed above, we find
%% it is close to the green(510nm)-white mixture line
%% we have constructed. 
selectInd = [];
for k = 1:size(munsell,2)
  [hueMod hueName lght sat] = parseMunsellName(munsellNames(k,:));
  if strcmp(hueName, 'G') & hueMod == 5
    selectInd = [ selectInd , k];
  end       
end
figure(2); hold on;
scatter(xyMunsImage(1,selectInd), xyMunsImage(2,selectInd), 5, 'b');
hold off;

% If you see a roughly parallel displacement between the '5G' munsell set
% and this green-white mixture line, it can be accounted for by the
% daylight used to illuminate the Munsell chips to form xyMunsImage.
% Using a 'white' illuminant with a constant spectral density
% instead, gives a better fit to the green-white mixture line.
XYZcoords = XYZ' * munsell;  % light sdf is constant, and is neglected
xyMunsCoords = XYZcoords(1:2,:) ./ repmat(sum(XYZcoords, 1), 2, 1);
xyMunsImage = xy2xyImage(xyMunsCoords);
figure(2); hold on;
scatter(xyMunsImage(1,selectInd), xyMunsImage(2,selectInd), 5, 'w');
hold off;
drawnow;

%% Close the CIE xy figure window
close(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COLOUR DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data for a Hitachi monitor.
load([iseRoot '/color/data/hit489.mat']);
gunSPD = measured;

% Spectral power distributions for the three phosphors
figure(1); clf;
plot(wavelength, gunSPD(:,3:-1:1));
drawnow;

% Any particular pixel on the display emits light with
% a SPD specified by a 3 vector of non-negative coefficients,
% called the gun coefficients, times the SPD for each of the
% three phosphors.  For example
coeffs = 0.33 * ones([3 1]);
displaySPD =  gunSPD * coeffs;
figure(1); hold on;
plot(wavelength, displaySPD, 'm');
hold off;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATING METAMERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppose we wish to simulate a particular SPD at a particular
% pixel.  For example, consider the Munsell reflectances
% illuminated by D65:
tDeg = 6.5;
daylight = daylightSpectrum(tDeg * 1000);
stimuli = munsell .* repmat(daylight, [1 size(munsell, 2)]);
% The first Munsell chip gives:
figure(1); clf;
plot(wavelength, stimuli(:,1));
drawnow;

% To get an perceptually indistinguishable source to one
% of these stimuli we need to ensure that the XYZ coeffs
% of the two stimuli are identical.  That is, we need
% to choose the gun coefficients, coeffs such that:
%%  XYZ' * displaySPD = XYZ' * stimuli
% But displaySPD is just phosphorSPD * coeffs, so:
%%  XYZ' * (gunSPD * coeffs) = XYZ'* stimuli
% Using the associativity of matrix multiplication:
%%  (XYZ' * gunSPD)  * coeffs = (XYZ' * stimuli)
% The first term on the left hand side is a 3 x 3 matrix
% formed from the XYZ coords of each of the phosphors:
hitXYZ = XYZ' * gunSPD
% The right hand side is the matrix of XYZ coefficients of each stimulus
stimXYZ = XYZ' * stimuli;
% For example the XYZ of the first stimulus is:
stimXYZ(:,1)
% The equation we need to solve for the coeffs is therefore:
%% hitXYZ * coeffs = stimXYZ

% Therefore we need to invert the hitXYZ matrix to determine Hitachi
% rgb gun coefficients in order to get the same XYZ value as the
% Munsell stimuli:
coeffs = hitXYZ \ stimXYZ;
% The first vector of gun coeffs is:
coeffs(:, 1)

% Sanity Check 1: 
% Compute the SPD of the light emitted by the display:
displaySPD = gunSPD * coeffs;
% Compare the XYZ coeffs of the emitted light from the display
% to the XYZ coeffs of the corresponding stimuli.  They
% should be the same:
[(XYZ' * displaySPD(:, 1))  stimXYZ(:, 1)]
% The maximum absolute error is:
fprintf(1, 'Sanity check 1: %e\n',...
            max(max(abs(stimXYZ - XYZ' * displaySPD))));
% which should be nearly zero.  

% The above sanity check shows that the
% selected coefficients cause the display to emit a light
% with the same XYZ coords, essentially, as the corresponding
% stimulus formed using a Munsell chip illuminated by D65.
% These stimuli will be perceptually indistinguishable by
% the standard observer.  

% Sanity Check 2:
% Compute the power absorbed by each of the cones for the
% first Munsell * D65 stimuli and for the first displaySPD:
xyzSensors' * [stimuli(:, 1) displaySPD(:,1)]
% These should be the same, indicating that they are perceptually
% indistinguishable.
% The maximum absolute error of the power captured by
% each of the cones, over all of the stimuli, is:
fprintf(1, 'Sanity check 2: %e\n',...
           max(max(abs(xyzSensors' * stimuli - xyzSensors' * displaySPD))));
% This should be nearly zero as well.


% However, the SPD for the display and the illuminated Munsell chip
% stimuli are NOT identical, even though their XYZ coefficients match.
% That is, they are metamers for each other.
figure(1); 
for k=selectInd
  clf;
  plot(wavelength, stimuli(:, k), 'g'); hold on;
  plot(wavelength, displaySPD(:, k), 'r');
  axis([mnW mxW 0 120]);
  title(['Metamer SPD for Munsell * D65 (g) and Display (r): Chip ' ...
      munsellNames(k,:)]);
  xlabel('Wavelength (nm)');
  ylabel('Power');
  hold off;
  pause(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COLOUR GAMUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A consequence of using three phosphors for the display is that
% the XYZ coefficients for any displayed pixel are linear
% combinations of the coefficients for each of the guns:
%% XYZ' * displaySPD = (XYZ' * gunSPD) * coeffs
% Moreover, the gun coefficients must all be non-negative.
% It then follows that the xy-coords for display pixels
% must lie within the triangle formed from the xy-coords for
% each of the phosphors:
xyGun = chromaticity( XYZ' * gunSPD);
xyImage = xy2xyImage(xyGun);

figure(1); clf;
[im map] = imread('data/diagxy.tif');
image(im); colormap(map); axis off;
title('CIE x,y Colour Space');
hold on;
for k=1:3
 j = mod(k, 3) + 1;
 x = [xyImage(1,k); xyImage(1,j)];
 y = [xyImage(2,k); xyImage(2,j)];
 line(x, y, 'Color', 'w');
end;
hold off;
drawnow;
% Thus the whole xy-space is not 'accessible' for display.
% Rather the range of colours that can be displayed, that is
% the colour gamut, is the interior of the triangle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Alternative Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([iseRoot '/color/data/sanyoLCD.mat']);
gunSPD = sanyoLCD;

% Superimpose the colour gamut for these basis spectral distributions.
xyGun = chromaticity( XYZ' * gunSPD);
xyImage = xy2xyImage(xyGun);
figure(1);
hold on;
for k=1:3
 j = mod(k, 3) + 1;
 x = [xyImage(1,k); xyImage(1,j)];
 y = [xyImage(2,k); xyImage(2,j)];
 line(x, y, 'Color', 'r');
end;
hold off;
drawnow;
% You can also redo the above calculations, starting just after the heading
% "COLOUR DISPLAY" and after the line in which gunSPD is set for
% the Hitachi display.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WHY USE JUST THREE CONE TYPES?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will show that the set of natural reflectances are well modelled
% by a three dimensional space.

% First set some preprocessing switches:  Leave them as
% they are set for the first pass through.  You may wish to
% change them later.  The results are not too sensitive on the
% settings used.
% for whether or not the process should
% subtracting out the mean of
% the data, using log reflectance, and/or a logistic model.
useMean=1;  % 1 => subtract out the mean reflectance
useLog=1;   % 1 => model the log reflectance instead of the reflectance.
            %    This avoids the need to worry about negative reflectances.
useLogistic = 0;  % First transform reflectance in (0,1)->(-infty, infty)
            % using the logistic map, r = refl/1.0-refl.  This ensures
            % that modelled reflectances will lie in the interval (0,1).

% We are going to fit a 3D linear model to a subset of the
% Munsell data, namely the data having indices in trainInd set below:
stimuli = munsell;
skip = 5;
trainInd = [1:skip:size(stimuli,2)];
size(trainInd,2)
% To use all the data, reset trainInd = [1:size(stimuli,2)]

if useLogistic
  stimuli = stimuli ./ (1.0-stimuli);
end
if useLog
     stimuli = log10(stimuli);
end

meanStimuli = sum(stimuli(:,trainInd), 2)/size(trainInd,2);
if useMean==1
  stimuli = stimuli - repmat(meanStimuli, [1 size(stimuli,2)]);
end

%% Plot the transformed reflectance distributions
figure(1); clf;
skip = 50;
skipInd = [1:skip:size(munsell,2)];
plot(wavelength, stimuli(:,skipInd));
title('Transformed Reflectances');
xlabel('Wavelength (nm)');
ylabel('Trans. Refl.');
hold on;
if (useMean ~=1)
  plot(wavelength, meanStimuli,'o');
end
hold off;
drawnow;

%% Compute the SVD of the training set (see svdTutorial.m)
trainStimuli = stimuli(:,trainInd);
[U, S, V] = svd(trainStimuli, 0);
sv = diag(S);
nSv = size(sv,1);

%% Plot the first 10 singular values
figure(1); clf;
k = [1:nSv];
kCrop = k(k<=10);
plot(kCrop, sv(k<=10)/max(sv), 'bo',...
     kCrop, sv(k<=10)/max(sv), 'b-');
axis([0 10 0 1]);
title('Relative size of singular values');
xlabel('Singular value index k');
ylabel('sv(k)/sv(1)');
drawnow;

% Compute the variance of each component, the total variance
% of the data set, and the cummulative variance of the first k
% components:
varComp = sv.^2 / size(trainInd,2);  
totalVar = sum(varComp);
q = cumsum(varComp)/totalVar;

%% Plot the proportion of the total variance accounted for
%% by the first k singular values:
rng = [1:10];
figure(1); clf;
plot(rng, q(rng>0), 'bo', rng, q(rng>0), 'b-');
title('Proportion of variance accounted for by k components');
xlabel('Singular value index k');
ylabel('Proportion of total variance');
drawnow;

% We see that just 3 linear components accounts for more than 95% of the
% total variance of the training data.  The precise fraction accounted
% for by 3 components is:
q(3)

% Another way to view the how the total variance of the training set
% is decomposed into the various components
% is to plot the percent of the total variance accounted
% for by each component separately:
figure(1); clf;
plot(rng, 100*varComp(rng>0)/totalVar, 'bo', ...
     rng, 100*varComp(rng>0)/totalVar, 'b-');
axis([0 10 0 20]);
title('Percent of variance accounted for by the k-th component');
xlabel('Singular value index k');
ylabel('Percent of total variance');
drawnow;

% Notice that the 3rd component accounts for more than 5% of the
% variance, but the 4th and higher components contribute less
% than 1 percent of the total variance each.
100 * varComp(rng>0)'/totalVar

%% The first four principal vectors encode, roughly speaking: 
%% 1: mean brightness,
%% 2: yellow vs blue, 
%% 3: green vs purple
%% 4,5,...: more rapid fluctuations of reflectance
figure(1); clf;
for k=1:10
  plot(wavelength, U(:,k));
  axis([mnW mxW -0.2 0.2]);
  hold on;
  line([mnW; mxW], [0;0], 'Color', 'k');
  title(sprintf('Principal reflectance component %d', k));
  xlabel('Wavelength (nm)');
  ylabel('Trans. Refl.');
  hold off;
  % fprintf(1, 'Type any key to continue...\n');
  pause(1);
end
% Notice that the higher principal vectors represent more
% rapid variations with respect to wavelength.  This is quite
% typical with principal components of natural stimuli.

%% The set of transformed Munsell reflectances looks like
%% an flat 3D elliptical pancake in a 361 dimensional space...
%% Show the scatter plot of the Munsell data projected
%% onto the first and cDim-th components, cDim = 2,3,...10.
%% There is 'internal structure' visible in the data set, especially
%% in the scatter plot for component 3 vs component 1. This
%% is due to the non-random nature of the Munsell set... after all
%% it was chosen to provide a systematic set of samples from colour space.
sigma = sqrt(varComp);
for cDim = 2:10
  figure(1); clf;
  plot(U(:,1)' * stimuli, U(:,cDim)' * stimuli, '.r'); axis equal;
  title(sprintf('Munsell data projected onto components 1 and %d', cDim));
  xlabel('Coeff. of 1st comp.');
  ylabel(sprintf('Coeff for comp. # %d', cDim));
  hold on;
  if useMean == 1 
    theta = [-1:0.01:1]*pi;
  else
    theta = [-1/2:0.01:1/2]*pi;
  end
  x = sigma(1)*cos(theta);
  y = sigma(cDim) * sin(theta);
  plot(x,y,'g', 2*x, 2*y, 'b');
  hold off;
  cDim
  % fprintf(1, 'Type any key to continue...\n');
  pause(1);
end
%% The curves show 1 and 2 standard deviation ellipses for
%% the Gaussian model, which provides a rough fit for the data.
%% These plots clearly demonstrate that the (transformed) Munsell
%% reflectances are approximated fairly well by a three dimensional
%% model.  The variation of the 4th and higher components is small
%% relative to the variation of the first component.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RECONSTRUCTING TRANSFORMED REFLECTANCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppose we use just the first three principal directions
nDim = 3;
pVec = U(:,[1:nDim]);

% The reconstructed reflectances are just the projections
% of the stimuli onto the subspace formed from the first three
% singular vectors:
recon = pVec * ( pVec' * stimuli);
skip = 50;
for k = 1:skip:size(stimuli,2)
  figure(1); clf; 
  plot(wavelength, stimuli(:,k), ':');
  hold on;
  plot(wavelength, recon(:,k), '-');
  hold off;
  title('Transformed refl (dotted), reconstruction(solid)');
  xlabel('Wavelength (nm)');
  ylabel('Trans. Refl.');
  pause(1);
end

% Compute some statistics for the errors in these reconstructions.
err = recon - stimuli;
sqrErr = sum(sum(err .* err));
totalVar = sum(sum(stimuli .* stimuli));
% Show the squared error in the reconstruction in terms of it's
% fraction of the total variance.  This should roughly agree
% with the 1 minus the fraction of the total variance accounted
% for by the first nDim components.  (The discrepancy will become
% smaller as the number of indices used for training (see trainInd
% above) is increased.)
[sqrErr/totalVar ,1-q(nDim)]

% Plot the distribution of errors in the reconstruction:
figure(1); clf;
[hist bin]=histo(reshape(err,prod(size(err)), 1), 64);
% Convert the histogram to a probability density
histDens = hist/prod(size(err))/(bin(2)-bin(1));
plot(bin, histDens);
title('Distribution of reconstruction errors: (data (b), fit (g) theory(r))');
xlabel('Error');
ylabel('Probability Density');
% Plot the Gaussian distribution with mean 0 and the sample error std dev 
hold on;
sigErr = std(reshape(err,prod(size(err)), 1))
dx = 0.1;
x = [-6:dx:6]*sigErr;
y = 1.0/sqrt(2.0*pi)/sigErr * exp(-(x .* x)/(2*sigErr*sigErr));
plot(x, y, 'g');
% Superimpose the Gaussian distribution with the std dev associated
% with the average variance of the error per pixel, based on the
% model using the singular values
k = [1:nSv];
varSigmaPixel = sum(varComp(k>nDim))/size(trainStimuli,1);
sigSigmaPixel = sqrt(varSigmaPixel)
dx = 0.1;
x = [-6:dx:6]*sigSigmaPixel;
y = 1.0/sqrt(2.0*pi)/sigSigmaPixel * ...
          exp(-(x .* x)/(2*sigSigmaPixel*sigSigmaPixel));
plot(x, y, 'r');
hold off;
% Notice the error distribution has a sharper peak near zero,
% and longer tails, than the two Gaussian models.
drawnow;


%% Histograms of coefficients are only roughly Gaussian, but
%% then again, the Munsell set is NOT a random sample.
tot = size(stimuli, 2);
for cDim=1:15
  [hist bin] = histo(U(:,cDim)' * stimuli/sigma(cDim), 64);
  figure(1); clf;
  histDens = hist/tot/(bin(2)-bin(1));
  % This normalizes the density so that:
  %    sum( histDens ) * (bin(2)-bin(1))
  % is one.
  plot(bin, histDens);
  axis([-6 6 0 0.8]);
  hold on;
  dx = 0.1;
  x = [-6:dx:6];
  y = 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2);
  plot(x, 1.0/sqrt(2.0*pi) * exp(-(x .* x)/2), 'g');
  % Here sum(y) * dx is one.
  title(sprintf('Coeff histogram for comp. # %d',cDim));
  pause(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RGB CORRELATIONS IN COLOUR IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We showed above that a large family of reflectances 
% are modelled very well with only a 3 principal components.
% Moreover, we showed that about 75% of the total variance
% was accounted for by the first component alone.

% A similar correlational structure should show up in
% colour images.  That is, we expect the brightness to provide
% the largest variation, with less variation, on average, in the
% yellow-blue components, and still less in the green-purple
% components.  Here we investigate whether or not this is
% the case.

% We have provided some colour images:
fn = {'hedvig3_050.tif'; % or use your own tif or tiff image.
      'allanA_080.tif';
      'hats.tiff'; 
      'lego-microserf.tif';
      'im7h.tif'};

% Loop over the image filenames in fn.
for kIm=1:size(fn,1)

  % Read the kIm^th image
  [im map] = imread(fn{kIm});

  % Convert the input image to a rgb image, where rgb(:,:,k)
  % refers to the red, green and blue channels for k=1,2,3, resp.
  % After conversion each rgb value is a double in the range [0, 1].
  if (size(size(im), 2) == 2)  % im is a 2D colour mapped image.
    sizeIm = size(im);
    rgbIm = map(double(im)+1, :);
    rgbIm = reshape(rgbIm, [sizeIm 3]);
  else
    rgbIm = double(im);
    if (max(max(rgbIm,1)) > 1.0)
      rgbIm = rgbIm/255;
    end
  end
  clear im map;

  examineRGB(rgbIm);  
  % This displays four figure windows.  Separate them on your screen.
  % See the comments in examineRGB.m for a description of their contents.

  fprintf(2, 'Press any key to continue...\n');
  pause;
end

