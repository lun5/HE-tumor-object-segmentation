function [hueMod, hueName, lght, satur] = parseMunsellName(name)
% function [hueMod, hueName, lght, satur] = parseMunsellName(name)
%
% Parse the name of a Munsell chip (eg. '5G 7/10') into:
%   - the hue is specified by hueMod and hueName:
%       hueMod either 2.5, 5, 7.5, 10
%       hueName either R, Y, G, B, or P or 
%               consecutive pairs YR, GY, ..., RP, where the order matters.
%   - the lightness, lght (a number),
%   - the saturation, satur (a number).

  [hueMod rem] = strtok(name,'RYGBP');
  hueMod = sscanf(hueMod, '%f');
  [hueName rem] = strtok(rem, ' ');
  [lght satur] = strtok(rem, '/');
  lght = sscanf(lght, '%f');
  satur = sscanf(satur(2:size(satur,2)), '%f');
  
