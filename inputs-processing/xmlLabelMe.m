%% modify the xml file from ImageScope to LabelMe
% Luong Nguyen
% 07/29/2015

%xmlfname = 'Z:\TilesForLabeling_bestImages\Annotation\2ALe5NgRyfnpo.xml';
%xmlfname = 'Z:\HEproject\data\Annotations\users\lun5\renamed_images\0anzqyibfuc.xml';
function new_xml = xmlLabelMe(xmlfname,folder_name,source)
[p,n,e] = fileparts(xmlfname);
if isempty(folder_name)
    folder_name = p;
end
[~, xml] = loadXML(xmlfname);
[~, tree] = xmlparse(xml);
new_xml = sprintf(['<annotation> \n <filename> %s.jpg </filename> \n <folder>' ...
    '%s</folder> \n <imagesize> \n <nrows>2048</nrows> \n' ...
    '<source>%s</source>\n'...
    '<ncols>2048</ncols>\n </imagesize>'],n,folder_name,source);

xml_lines = tree.tagname;
IndexC = strfind(xml_lines,'Region Id=');
Index = find(not(cellfun('isempty', IndexC)));
IndexC = strfind(xml_lines,'Plots');
End_Index = find(not(cellfun('isempty',IndexC)));
for i = 1:length(Index)
    s = strsplit(xml_lines{Index(i)},' ');
    indx = strfind(s,'Text=');
    indx = find(not(cellfun('isempty',indx)));
    object_name = strsplit(strrep(s{indx},'"',''),'=');
    object_name = strcat(object_name{2:end});
    xml_str = sprintf(['\n<object>\n <name>%s</name>\n<deleted>0</deleted>\n'...
        '<verified>0</verified>\n <occluded>no</occluded>\n <attributes />\n' ...
        '<parts><hasparts /> <ispartof /> </parts>\n<date></date>\n'...
        '<id>%d</id>\n<polygon>\n<username>%s</username>\n'],object_name,i-1,source);
    if i < length(Index)
        end_indx = Index(i+1)-1;
    else
        end_indx = End_Index -1;
    end
    for j = (Index(i)+3) : end_indx
       coords = xml_lines{j};
       coords = strrep(coords,'"','');
       coords = strrep(coords,'Vertex ','');
       coords = strsplit(coords,' ');
       xcoord = strsplit(coords{1},'=');
       xcoord = str2double(xcoord{2});
       ycoord = strsplit(coords{2},'=');
       ycoord = str2double(ycoord{2});
       coord_str = sprintf('\n<pt>\n<x>%d</x>\n<y>%d</y>\n</pt>\n',round(xcoord),round(ycoord));
       xml_str = strcat(xml_str,coord_str);
    end
    new_xml = strcat(new_xml,xml_str,sprintf('\n</polygon>\n</object>\n'));
end

new_xml = strcat(new_xml,sprintf('\n</annotation>'));
%[v,tree] = xmlparse(new_xml);
end