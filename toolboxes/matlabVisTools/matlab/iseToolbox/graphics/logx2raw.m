function logx2raw(base,precision)
%logx2raw([base],[precision])
%
%Converts X-axis labels from log to raw values.
%base:    	base of log transform; default base is e.
%precision:	number of decimal places;  default  is two.
%
%Example:
% x=linspace(-3,0,11);
% plot(log(x),log(x.^2));
% logx2raw
% logy2raw

%SEE ALSO;   Logy2raw
%11/17/96	gmb	Wrote it.
%6/6/96	        gmb added precision argument

if nargin==0
	base=exp(1);
end
if nargin<2
	precision=2;
end
qt='''';
loglab=get(gca,'XTickLabel');
maxlen=0;
for i=1:size(loglab,1);
	rawct=base.^(str2num(loglab(i,:)));
	estr=(['xtli=num2str(sprintf(',qt,'%',sprintf('%2.1f',precision*1.1), ...
	'f',qt,',rawct));']);
	eval(estr);
	%xtli=num2str(sprintf('%3.3f',rawct));
	if (length(xtli)>maxlen)
		maxlen=length(xtli);
	end
	eval(['xtl',num2str(i),' = xtli;'])
end

for i=1:size(loglab,1);
	eval(['xtli=xtl',num2str(i),';']);
	for j=1:maxlen-size(xtli,2)
		xtl(i,j)=' ';
	end
	if (size(xtli,2)>0)
		xtl(i,maxlen-length(xtli)+1:maxlen)=xtli;
	end
end
set(gca,'XTickLabel',xtl);
