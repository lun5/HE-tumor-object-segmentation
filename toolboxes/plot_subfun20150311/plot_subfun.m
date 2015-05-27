function varargout = plot_subfun(foo,varargin)
%plot_subfun(foo)
%plots a dependency tree of the subroutines within a function file
%
%EXAMPLES
%plot_subfun(foo,'-hide','a','b','c')
%hides subroutines 'a','b','c'
%also hides their "kids" (any functions that are only called by a,b,c)
%
%plot_subfun(foo,'-kids','a','b',c')
%as above, but only hides the kids, not the functions themselves
%
%the flags '-hide' and '-kids' affect the names that follow
%the default is '-hide', so:
%plot_subfun('a','b','c','-kids','d','e','f','-hide','g','h')
%will apply 'hide' to a,b,c,g,h and 'kids' to d,e,f
%
%it is possible to specify the names of functions that are not present
%these will simply be ignored
%
%ALL HIDE OPTIONS
%-unused        show unused functions too
%-unusedonly    show only unused functions
%-hide          hide these functions and their kids
%-kids          hide the kids of these functions
%-nest          hide nested functions   (those within another function)
%-ends          hide end functions      (have one parent and no kids)
%
%these options are applied in the order shown above. functions that
%have been hidden by a previous option are treated as not existing for the
%purpose of subsequent options.
%
%OTHER OPTIONS
%-trim          trim hidden functions from output
%-list          show list of subroutines in command window
%-extlist       show list of external functions called
%-extshow       external function calls included on figure
%-ext           ditto
%-noplot        do not plot the figure
%
%The author wishes to acknowledge the following file exchange submissions:
%farg, ftoc, fdep, mGraphViz, mkdotfile

data = sub_setup(foo,varargin{:});  %parse input and find dependencies
data = sub_hide(data);              %hide functions
data.fig = sub_show_plot(data);     %plot the figure
sub_show_list(data);                %show a list of functions
if nargout;                         %return output only if requested
    data = sub_trim(data);          %trim out hidden functions
    varargout{1} = data;
end
end

%% SETUP
function data = sub_setup(foo,varargin)
sub_check_options(varargin{:})
data = sub_parse_options(varargin{:});  %parse options
data = sub_parse_file(data,foo);        %parse function structure
data = sub_extl_calls(data);            %remove external calls
data = sub_deps(data);                  %find dependencies from remaining calls
data = sub_find_used(data);             %find functions that are actually used
data = sub_unique_names(data);          %change function names to all be unique

function sub_check_options(varargin)
%check input options

%check all options are strings or cellstring arrays
if ~all(cellfun(@ischar,varargin) | cellfun(@iscellstr,varargin))
    error('all inputs must be strings or cellstring arrays');
end

%any starting with - must be a known keyword
str = varargin(~cellfun('isempty',regexp(varargin,'^-.*')));
targ = {...
    '-unused','-unusedonly','-hide','-kids','-nest','-ends',...
    '-trim','-list','-ext','-noplot','-unhide','-extlist','-extshow'};
tf = ismember(str,targ);
if ~all(tf);
    disp(str(~tf));
    warning('the above are not recognised options, they will be ignored');
end

end

function data = sub_parse_options(varargin)

%parse options
data.out.extlist= false;    %show external functions as list
data.out.extshow= false;    %show external functions on figure
data.out.list   = false;    %display list
data.out.plot   = true;     %display figure
data.out.trim   = false;    %trim hidden functions from output
data.hide.nest  = false;    %hide nested functions
data.hide.ends  = false;    %hide deadend functions
data.hide.orph  = true;     %hide orphan functions
data.hide.used  = false;
data.hide.funs  = {};       %functions to hide with kids
data.hide.kids  = {};       %functions to hide kids only
data.out.unhide  = false;    %show functions that have been hidden
addto = 'funs'; %file names should be added to ignore or nokid
for i=1:numel(varargin)
    val = varargin{i};
    switch val
        case '-hide';       addto = 'funs';
        case '-kids';       addto = 'kids';
        case '-ext' ;       data.out.extshow= true;
        case '-extshow';    data.out.extshow= true;
        case '-extlist';    data.out.extlist= true;
        case '-list';       data.out.list   = true;
        case '-noplot';     data.out.plot   = false;
        case '-trim';       data.out.trim   = true;
        case '-nest';       data.hide.nest  = true;
        case '-ends';       data.hide.ends  = true;
        case '-unused';     data.hide.orph  = false;
        case '-unhide';     data.out.unhide = true;
        case '-unusedonly'; data.hide.orph = false; 
                            data.hide.used = true;
        otherwise
            data.hide.(addto) = [data.hide.(addto) val];
    end
end

%hide automatically includes kids
data.hide.funs = unique(data.hide.funs);
data.hide.kids = unique([data.hide.kids data.hide.funs]);
end

function data = sub_parse_file(data,foo)
%process the file
data.foo = which(foo);
%20150302 catch built-in function
if regexp(data.foo,'^built-in')
    data.foo = regexprep(data.foo,'^built-in\s+\(','');
    data.foo = regexprep(data.foo,'\)$','');
    if isempty(regexp(data.foo,'\.m$','once')); data.foo = [data.foo '.m']; end
end

if ~exist(data.foo,'file');
    error('plot_subfun:noFile',['File "' data.foo '" does not exist']);
end
p=mlintmex(data.foo,'-m3','-calls');
b = regexp(p,'(?<kind>(M|S|N))(?<gen>\d) (?<line>\d+) \d+ (?<name>\w+)','names');
%end of functions
e = regexp(p,'E. (?<line>\d+) \d+ (?<name>\w+)','names');
%all calls in functions
c = regexp(p,'U. (?<line>\d+) \d+ (?<name>\w+)','names');

%2015.01.30 create an empty set of functions as a placeholder
%this should cope with scripts that contain no subroutines
data.fun.name = {};
data.fun.beg = [];
data.fun.end = [];
data.fun.gen = [];
data.fun.kind = '';

for i=1:numel(b)
    data.fun.name{i} = b(i).name;
    data.fun.beg(i) = str2double(b(i).line);
    data.fun.end(i) = str2double(e(i).line);
    data.fun.gen(i) = str2double(b(i).gen);
    data.fun.kind(i) = b(i).kind;
end
for i=1:numel(c)
    data.call.name{i} = c(i).name;
    data.call.line(i) = str2double(c(i).line);
end

%catch empty cases
if numel(b);
    data.fun.numof = numel(data.fun.name);
    data.fun.show = true(1,data.fun.numof);
    data.fun.deps = false(data.fun.numof,data.fun.numof);
    data.fun.used = true(1,data.fun.numof);
    data.fun.hide = false(1,data.fun.numof);
else
    data.fun.numof = 0;
    data.fun.show = [];
    data.fun.deps = [];
    data.fun.used = [];
    data.fun.hide = [];
end

if ~numel(c);
    data.call.name = {};
    data.call.line = [];
end


end

function data = sub_extl_calls(data)
%data = sub_extl_calls(data);

%% handle external calls

%which calls are to external functions
data.call.isext = ~ismember(data.call.name,data.fun.name);

%make list of external calls
names = unique(data.call.name(data.call.isext));
%remove calls to built-in functions
a = cellfun(@which,names,'UniformOutput',false);
a = regexp(strrep(a,'\','/'),strrep(matlabroot,'\','/'),'match','once');
keep = cellfun(@isempty,a);
data.external = names(keep);

%If external functions not shown on the figure, ignore these calls
if ~data.out.extshow
    %remove from list of calls
    data.call.name  = data.call.name( ~data.call.isext);
    data.call.line  = data.call.line( ~data.call.isext);
    data.call.isext = data.call.isext(~data.call.isext);
end

%if external functions shown on the figure, must add them to the list of
%functions (with -1 line number, or additional field to flag them as
%external ?)
if data.out.extshow
    ne = numel(data.external);
    nf = data.fun.numof;
    data.fun.name = [data.fun.name data.external];
    data.fun.beg  = [data.fun.beg  nan(1,ne)];
    data.fun.end  = [data.fun.end  nan(1,ne)];
    data.fun.gen  = [data.fun.gen  zeros(1,ne)];
    data.fun.kind = [data.fun.kind repmat('E',[1 ne])];
    data.fun.numof = data.fun.numof + ne;
    data.fun.show = [data.fun.show ones(1,ne)];
    data.fun.deps = zeros(ne+nf,ne+nf);
    data.fun.used = [data.fun.used true(1,ne)];
    data.fun.hide = [data.fun.hide false(1,ne)];
end

end

function data = sub_deps(data)
for i=1:numel(data.fun.name)
    %lines is lines owned by this function only
    lines = data.fun.beg(i):data.fun.end(i);
    %remove any lines from nested functions within this function
    %all lines belonging to higher generation functions
    for I = find(data.fun.gen>data.fun.gen(i));
        lines(lines>=data.fun.beg(I) & lines<=data.fun.end(I)) = [];
    end
    called = unique(data.call.name(ismember(data.call.line,lines)));
    
    data.fun.deps(i,:) = false(size(data.fun.name));
    for j=1:numel(called)
        tf = ismember(data.fun.name,called{j});
        if sum(tf)==1;
            data.fun.deps(i,tf) = true;
            continue;
        end
        if sum(tf)>1;
            %need to find the one of several options being called
            I = find(tf);
            
            %recursively search up through my parents, starting with me
            done = false;
            me = i;
            while ~done;
                %ignore any that are above gen(me)+1
                I(data.fun.gen(I)>data.fun.gen(me)+1) = [];
                %if it is in me and my gen+1, use this
                isin = ...
                    data.fun.beg(I)>data.fun.beg(me) & ...
                    data.fun.end(I)<data.fun.end(me);
                ischild = data.fun.gen(I) == data.fun.gen(me)+1 & isin;
                if any(ischild);
                    data.fun.deps(i,I(ischild)) = true;
                    done = true;
                    continue
                end
                
                %if we are gen=0, are there any other gen=0 that match ?
                if data.fun.gen(me)==0 && any(data.fun.gen(I) == 0);
                    I(data.fun.gen(I)>0) = [];
                    data.fun.deps(i,I) = true;
                    done = true;
                    continue;
                end
                
                %recursively search up through my parents doing the same thing
                if data.fun.gen(me)>0; me = sub_dad(data,me); continue; end
                
                %if we are here, we are gen 0 and haven't found it, give up
                keyboard
                done = true;
            end
        end
    end
end
%TODO deal with non-unique funciton names
%TODO list plot : plot text to figure, draw lines on that.

    function dad = sub_dad(data,me)
        dad = find(...
            data.fun.beg<data.fun.beg(me) & ...
            data.fun.end>data.fun.end(me) & ...
            data.fun.gen==data.fun.gen(me)-1);
        if isempty(dad) || numel(dad)>1;
            keyboard
        end
    end
end

function data = sub_find_used(data)
%flag functions as "used" if main function depends on them

%catch functions that are empty
if ~data.fun.numof; return; end

%for the sake of these checks, ignore calls to self
deps = data.fun.deps; deps(eye(size(deps))==1) = false;

%for the sake of these checks, hidden functions have no kids/parents
deps(~data.fun.show,:) = false;
deps(:,~data.fun.show) = false;

%start with only first (main) function shown
used = false(size(data.fun.name));
used(1) = true;

%recursively : all children of shown, add to set of shown
used_old = 0;
while ~isequal(used,used_old)
    used_old = used;
    used = used | any(deps(used,:),1);
end

%update functions shown
data.fun.used = used;

end

function data = sub_unique_names(data)
%invents unique names (preferably of the same length as the original) for
%each function. This is done so different subroutines with the same name
%are shown correctly, rather than being lumped together as one function

%20150130 catch empty
if isempty(data.fun.name); return; end

%TODO20150215 tag names changed below, so they get changed back!
%do not rely on string-replace, as that may change functionnames that
%really were sub_a__line231_
for names = unique(data.fun.name)
    name = names{1};
    I = find(ismember(data.fun.name,name));
    if numel(I)>1
        for i = I
            data.fun.name{i} = [data.fun.name{i} '__line' num2str(data.fun.beg(i)) '_'];
        end
    end
end
end
end

%% HIDE FUNCTIONS

function data = sub_hide(data)
data = sub_hide_used(data);         %decide if used/unused are shown
data = sub_hide_funs(data);         %hide functions
data = sub_hide_kids(data);         %hide kids
data = sub_hide_nest(data);         %hide nested
data = sub_hide_ends(data);         %hide one-parent dead ends
%data = sub_hide_orphan(data);       %hide orphans
end

function data = sub_hide_used(data)
show = true(size(data.fun.name));
if data.hide.used; show( data.fun.used) = false; end
if data.hide.orph; show(~data.fun.used) = false; end
data.fun.show = data.fun.show & show;
end

function data = sub_hide_funs(data)
if isempty(data.hide.funs); return; end

[tf,loc] = ismember(data.hide.funs,data.fun.name);
loc = loc(tf);
show = true(size(data.fun.name));
show(loc) = false;
data.fun.show = data.fun.show & show;
data.fun.hide(~show) = true;
end

function data = sub_hide_kids(data)
if isempty(data.hide.kids); return; end

%local version of dependencies, we can muck with. do not feed up to parent.
%for the purpose of these checks, links to self do not count
deps = data.fun.deps;
deps(eye(size(deps))==1) = false;

%names of functions whose kids to hide
names = unique([data.hide.kids data.hide.funs]);

%functions that originally had no dads, these will not be hidden
orig = sum(deps,1)==0;

%also, will never hide themselves when doing the kids thing
orig(ismember(data.fun.name,names)) = true;

nokids_old = nan;
nokids_new = ismember(data.fun.name,names);
%iterate until nokids does not change
while ~isequal(nokids_old,nokids_new)
    %sever connections to their kids
    deps(nokids_new,:) = false;
    %update set nodads (ignore orig)
    nodads = sum(deps,1)==0;
    nodads(orig) = false;
    %these are now nokids
    nokids_old = nokids_new;
    nokids_new = nodads;
end
data.fun.show(nokids_new) = false;
data.fun.hide(nokids_new) = true;
end

function data = sub_hide_nest(data)
if ~data.hide.nest; return; end %if nested option not received, do nothing
nest = data.fun.gen>0; %functions that are nested
for i=find(nest)
    if ~data.fun.show(i); continue; end %if already hidden, ignore
    %my parents inherit my children
    dads = data.fun.deps(:,i);
    data.fun.deps(dads,:) = data.fun.deps(dads,:) | repmat(data.fun.deps(i,:),[sum(dads) 1]);
    %ignore me
    data.fun.show(i) = false;
    data.fun.hide(i) = true;
end
end

function data = sub_hide_ends(data)
if ~data.hide.ends; return; end

%for the purpose of this test, calls to self to not count
deps = data.fun.deps;
deps(eye(size(deps))==1) = false;

%for the purpose of this test, hidden functions do not count
deps(~data.fun.show,:) = false;
deps(:,~data.fun.show) = false;

%hide dead ends, with one parent and no kids
nokids = ~any(deps,2)';
onedad = sum(deps,1)==1;
data.fun.show(nokids & onedad) = false;
data.fun.hide(nokids & onedad) = true;
end

%% PLOT, DISPLAY AND OUTPUT
function opt = sub_show_plot(data)
data = sub_trim(data,true);
data = sub_add_fake_self(data);

if ~data.out.plot; opt = []; return; end
if ~any(data.fun.deps(:)); disp('There are no subroutines to show'); opt=[]; return; end

%colours of each box
deps = data.fun.deps;
deps(eye(size(deps))==1) = false; %for the purpose of this check (only) ignore calls to self
I = ~any(deps,1); %has no parents
J = data.fun.gen>0; %is nested
K = data.fun.kind == 'E';
cols = repmat('b',size(data.fun.name));
cols(K) = 'k';
cols(J) = 'g';
cols(I) = 'r';

%shape of each box (rounded or square corners)
shaps(~K) = 'r';
shaps(K)  = 's';

%convert the dependencies into from/to pairs
[from,to] = ind2sub(size(deps),find(data.fun.deps));
opt = plot_graph(data.fun.name,from,to,'-colours',cols,'-shape',shaps);

%if external called, add a list of external function to the right hand side
%of the plot
if data.out.extlist
    %add 20% to the right axes limits
    x = get(gca,'xLim');
    y = get(gca,'yLim');
    dx = [0 0.2*diff(x)];
    set(gca,'xLim',x+dx);
    %write into this space
    txt = data.external; if isempty(txt); txt = {'none'}; end
    txt = [{'EXTERNAL CALLS:','=================='},txt];
    txt = strrep(txt,'_','\_');
    %text(x(2),y(2)-0.1*diff(y),txt,...
    %    'VerticalAlignment','Top');
    text(x(2),y(2),txt,...
        'VerticalAlignment','Top');
end

%correct names of type "sub_a__line203_" to "sub_a (line203)"
opt = sub_show_fix_names(opt);

[opt,data] = sub_remove_fake_self(opt,data);
opt = sub_show_unhide(opt,data); %change colour of -unhide function to grey
end

function opt = sub_show_unhide(opt,data)
%boxes forcibly shown with the -unhide option are coloured grey as it their
%text and any lines to/from them

%which function are hidden
if ~data.out.unhide; return; end
targ = find(data.fun.hide);

%find hidden nodes, colour them grey
tf = ismember([opt.node(:).index],targ); %nodes that are hidden
for i=find(tf)
    set(opt.node(i).handle_box,'EdgeColor',[.9 .9 .9])
    set(opt.node(i).handle_text,'Color',[.9 .9 .9]);
end

%edges to/from hidden boxes also coloured grey
val = vertcat(opt.edge.index);
tf = ismember(val(:,1),targ) | ismember(val(:,2),targ);
for i=find(tf')
    set(opt.edge(i).handle,'Color',[.9 .9 .9]);
    set(opt.edge(i).handle_arrow,'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9]);
end
end

function opt = sub_show_fix_names(opt)
%correct names of type "sub_a__line203_" to "sub_a (line203)"
for i=1:numel(opt.node)
    str = opt.node(i).text;
    if regexp(str,'__line\d+_$');
        str = regexprep(str,'__line',' (line');
        str = regexprep(str,'_$',')');
        set(opt.node(i).handle_text,'String',str);
    end
end
end

function data = sub_add_fake_self(data)
deps = data.fun.deps;
%index of the self references
I = eye(size(deps))==1;
%the ones that don't have a self reference become fake
data.fakeself = find(~deps(I));
%create these as links to be removed later
I = sub2ind(size(deps),data.fakeself,data.fakeself);
data.fun.deps(I) = true;
end

function [opt,data] = sub_remove_fake_self(opt,data)

%find list of function names that have fake selfs:
names = data.fun.name(data.fakeself);
keep = true(size(opt.edge));
for i=1:numel(opt.edge)
    e = opt.edge(i);
    if ~isequal(e.endpoints{1},e.endpoints{2}); continue; end
    if ismember(e.endpoints{1},names);
        keep(i) = false;
    end
end

%remove the fake edges from the figure
for I = find(~keep);
    delete(opt.edge(I).handle);
    delete(opt.edge(I).handle_arrow);
end
%remove the fake edges from the figure
opt.edge = opt.edge(keep);
data = rmfield(data,'fakeself');
end

function sub_show_list(data)
if ~data.out.list; return; end

for i=1:numel(data.fun.name)
    name = data.fun.name{i};
    ibeg = data.fun.beg(i);
    iend = data.fun.end(i);
    fprintf(' \n');
    fprintf('%s (line %d-%d)\n',name,ibeg,iend);
    for j=find(data.fun.deps(i,:))
        fprintf('-%s\n',data.fun.name{j})
    end
end
end

function data = sub_trim(data,varargin)
if ~data.out.trim && nargin<2; return; end
if data.out.unhide
    keep = data.fun.show | data.fun.hide;
else
    keep = data.fun.show;
end

for name = {'name','beg','end','gen','kind','show','used','hide'}
    val = data.fun.(name{1});
    data.fun.(name{1}) = val(keep);
end
data.fun.deps = data.fun.deps(keep,:);
data.fun.deps = data.fun.deps(:,keep);

data.fun.numof = sum(keep);

%remove calls to removed functions
keep = ismember(data.call.name,data.fun.name);
data.call.name = data.call.name(keep);
data.call.line = data.call.line(keep);
end

%% SUBROUTINES THAT DO NOTHING
%These are only here to be seen in the codetree of this function
%to show how looped and orphaned functions look.

function sub_a() %#ok<DEFNU>
sub_aa()
    function sub_aa()
        sub_b()
    end
end

function sub_b()
sub_a()
    function sub_a()
        sub_b()
        function sub_b()
        end
    end
end

function sub_notcalled() %#ok<DEFNU>
sub_selfcalled()
sub_nested()
    function sub_nested()
    end
end

function sub_selfcalled()
sub_selfcalled()
end

function sub_pair1()
sub_pair2()
end

function sub_pair2()
sub_pair1()
sub_deadend()
end

function sub_deadend()
end

function twodad1()
twodad3()
end

function twodad2() %#ok<DEFNU>
twodad3()
end

function twodad3()
twodad4()
end

function twodad4()
twodad5;
end

function twodad5()
twodad1;
end

function isolatedloop()
isolatedloop();
end

%% DEVNOTES
%20140901   added handling of multiple instances of the same function name
%20150215   cleaned up display of multiple subrotines that share a name:
%           "sub_a__line231_" becomes "sub_a (line231)"
%DONE       handles subroutines truly called sub_a__line231_
%DONE       add '-kids' option, that removes children instead of self.
%DONE       if you want to show a function with no deps : make deps to self, then remove link on plot
%20150228   rewritten based on user feedback
%           reorganised code for modularity
%           added additional options for show/hide
%20150301   bugfix to handle built-in function
%           bugfix to handle empty functions
%DONE       -unhide option for demos, shows in grey (smaller ?)
%TODO       sizing by number of lines, number of dependents, user input
%TODO       documentations
%TODO       colorblind mode (custom colours)
%TODO       warn when analysing a built-in function
%TODO       deal with callback creation (Joshua was reporting a problem)
%TODO       extlist vs extshow options (ext also -> extlist)
%           ext to show as black. list separate
%TODO       builtin option ?
