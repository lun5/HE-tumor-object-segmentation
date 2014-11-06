function scaleYaxis(k,precision)
%scaleYaxis(k,precision)
%
%Scales Y-axis labels by a factor of k
%precision:	number of decimal places;  default  is two.
%
%SEE ALSO;   Logy2raw

%1/10/97    gmb wrote scaleYaxis from logy2raw.

if nargin==0
	k=1;
end
if nargin<2
	precision=2;
end
qt='''';
lab=get(gca,'YTickLabels');
maxlen=0;
for i=1:size(lab,1);
	rawct=k*(str2num(lab(i,:)));
	rawct=k*(str2num(lab(i,:)));
	estr=(['xtli=num2str(sprintf(',qt,'%',sprintf('%2.1f',precision*1.1), ...
	'f',qt,',rawct));']);
	eval(estr);
	%xtli=num2str(sprintf('%3.3f',rawct));
	if (length(xtli)>maxlen)
		maxlen=length(xtli);
	end
	eval(['xtl',num2str(i),' = xtli;'])
%	xtl(i,:)=num2str(sprintf('%2.2f',rawct));
end

for i=1:size(lab,1);
	eval(['xtli=xtl',num2str(i),';']);
	for j=1:maxlen-size(xtli,2)
		xtl(i,j)=' ';
	end
	if (size(xtli,2)>0)
		xtl(i,maxlen-length(xtli)+1:maxlen)=xtli;
	end
end
set(gca,'YTickLabels',xtl);


