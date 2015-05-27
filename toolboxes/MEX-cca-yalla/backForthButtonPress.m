function [k, stop, varargout] = backForthButtonPress(h, k, minK, maxK, cnt) 
%%%% [k, stop, [cnt]] = backForthButtonPress(h, k, minK, maxK, cnt) 
     
% Increment/decrement index k based on the key/button press in 
% the figure window. 
%   h = figure handle 
%   k = current index 
%   minK <= k <= maxK  range of index. 
%   cnt  number of calls to backForthButton press in which 
%        to print out the user prompt. 
% Keyboard Controls: 
%   if ENTER > . SPACEBAR or LeftMouseButton is pressed 
%      k increases by 1, or sticks at maxK  
%   if BKSPC < or , is pressed 
%      k decreases by 1, or sticks at minK  
%   if CTRL_ENTER is pressed 
%      k reset to maxK 
%   if CTRL_BKSPC is pressed 
%      k reset to minK 
%   if q, Q, ESC is pressed 
%      old k, with stop = 1 is returned 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  if nargin<5 
    cnt = 0; 
  end 
 
  PROMPT = ['\nMOVE THE CURSOR ONTO THE FIGURE WINDOW AND TYPE:\n' ... 
            ' Enter Bkspc Cntrl_Enter Cntrl_Bkspc ESC  ?\n' ... 
            ' next  prev  first       last        quit this\n']; 
  ENTER = char(13); 
  ESC = char(27); 
  SPACEBAR = char(32); 
  BKSPC = char(8); 
  CTRL_ENTER = char(10); 
  CTRL_BKSPC = char(127); 
  if cnt > 0 
    fprintf(2,PROMPT); 
    cnt = cnt-1; 
  end 
  varargout(1) = {cnt}; 
 
  if minK > maxK 
    error('BackForthButtonPress: requires minK <= maxK'); 
  end 
  stop = 0; 
  if 0 == waitforbuttonpress  % Mouse button press 
    k = max(k + 1, minK); 
    k = min(k, maxK); 
  else     % Keyboard button press 
    c = get(h,'CurrentCharacter');
    %double(c) % Dump character as an int 
    if c == ENTER | c == '.' | c == SPACEBAR 
      k = min(k+1, maxK); 
    elseif c == CTRL_ENTER  | c == '>' 
      k = maxK; 
    elseif c == BKSPC | c == ','   
      k = max(k-1, minK); 
    elseif c == CTRL_BKSPC  | c == '<'  
      k = minK; 
    elseif  c == 'q' | c == 'Q' | c == ESC 
      stop = 1; 
    elseif c == '?' 
      fprintf(2, PROMPT); 
    end 
  end 
