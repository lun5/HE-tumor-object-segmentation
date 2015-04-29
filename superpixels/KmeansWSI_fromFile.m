

function KmeansWSI_fromFile(fullFileName,k)


indd=strfind(fullFileName, '/');
indx=indd(end);
imageInfo = imfinfo(fullFileName);
firstLayer_info = imageInfo(1,1);



ww = firstLayer_info.Width;
ww = round(ww/4000);
hh = firstLayer_info.Height;
hh = round(hh/4000);
delw = round(firstLayer_info.Width/ww);
delh = round(firstLayer_info.Height/hh);

temp_WSI = imread(fullFileName);

scaled_image= imresize(temp_WSI,0.25);
nrows  = size(scaled_image,1);
ncols= size(scaled_image,2);
%RGB
Data = double([reshape(scaled_image(:,:,1),nrows*ncols,1) reshape(scaled_image(:,:,2),nrows*ncols,1) reshape(scaled_image(:,:,3),nrows*ncols,1)]);
%HSV
% rotation_matrix=[0.595411526389924,0.535999946816847,0.598489073629885;-0.349245071618676,-0.498189188484696,0.793621706121390;-0.723541978182390,0.681550870800777,0.109432245332974];
% r = in_image(:,:,1)./255; g = in_image(:,:,2)./255; b = in_image(:,:,3)./255;
% rotated_coordinates = rotation_matrix*double([r(:)'; g(:)'; b(:)']);
% theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
% Data = theta';


kMeansClusterCenters = initializeClusterVectorsUsingPrincipalComponent(Data,k);
m= kMeansClusterCenters;
% determine lumen,strome, nuclei labels RGB
if ((m(1, 1) > m(2, 1)) && (m(1, 1) > m(3, 1)))
    lm_id = 1;
    if (m(2, 1) > m(3, 1))
        st_id = 2;
        nc_id = 3;
    else
        st_id = 3;
        nc_id = 2;
    end
else
    if ((m(2, 1) > m(1, 1)) && (m(2, 1) > m(3, 1)))
        lm_id = 2;
        if (m(1, 1) > m(3, 1))
            st_id = 1;
            nc_id = 3;
        else
            st_id = 3;
            nc_id = 1;
        end
    else
        lm_id = 3;
        if (m(1, 1) > m(2, 1))
            st_id = 1;
            nc_id = 2;
        else
            st_id = 2;
            nc_id = 1;
        end
    end
end

% %determine lumen,strome, nuclei labels Hue
% if ((m(1, 1) > m(2, 1)) && (m(1, 1) > m(3, 1)))
%     lm_id = 1;
%     if (m(2, 1) > m(3, 1))
%         st_id = 3;
%         nc_id = 2;
%     else
%         st_id = 2;
%         nc_id = 3;
%     end
% else
%     if ((m(2, 1) > m(1, 1)) && (m(2, 1) > m(3, 1)))
%         lm_id = 2;
%         if (m(1, 1) > m(3, 1))
%             st_id = 3;
%             nc_id = 1;
%         else
%             st_id = 1;
%             nc_id = 3;
%         end
%     else
%         lm_id = 3;
%         if (m(1, 1) > m(2, 1))
%             st_id = 2;
%             nc_id = 1;
%         else
%             st_id = 1;
%             nc_id = 2;
%         end
%     end
% end

testIndX=1;
testIndY=1;
for ii=1:delw:firstLayer_info.Width
    for jj=1:delh:firstLayer_info.Height
        endi=ii+delw;
        endj=jj+delh;
        if ii+delw>firstLayer_info.Width
            endi=firstLayer_info.Width;
        end
        if jj+delh>firstLayer_info.Height
            endj=firstLayer_info.Height;
        end
        if abs(endi-ii)<50
            continue;
        end
        if abs(endj-jj)<50
            continue;
        end
        KmeansImFileName = [fullFileName(indx+1:end-4) '_ImX' num2str(testIndX) '_ImY' num2str(testIndY) '_stX' num2str(ii) '_stY' num2str(jj) '_k3'];
        
        if exist(KmeansImFileName,'file')>0
          testIndY=testIndY+1;
          continue;
        end

        tmpim= imread([pwd '/' fullFileName(indx+1:end)],'PixelRegion', {[jj endj] , [ii endi]});
        imwrite(tmpim,[fullFileName(indx+1:end-4) '_ImX' num2str(testIndX) '_ImY' num2str(testIndY) '_stX' num2str(ii) '_stY' num2str(jj) '.jpg'],'jpg');
        
        mrows  = size(tmpim,1);
        mcols= size(tmpim,2);
        %RGB
                Datam = double([reshape(tmpim(:,:,1),mrows*mcols,1) reshape(tmpim(:,:,2),mrows*mcols,1) reshape(tmpim(:,:,3),mrows*mcols,1)]);
        %HSV
%         r = tmpim(:,:,1)./255; g = tmpim(:,:,2)./255; b = tmpim(:,:,3)./255;
%         rotated_coordinates = rotation_matrix*double([r(:)'; g(:)'; b(:)']);
%         theta = angle(rotated_coordinates(2,:) + 1i*rotated_coordinates(3,:));
%         
%         im_theta = reshape(theta,size(r));
%         imwrite(im_theta,[fullFileName(indx+1:end-4) '_ImX' num2str(testIndX) '_ImY' num2str(testIndY) '_stX' num2str(ii) '_stY' num2str(jj) '.jpg'],'jpg');
%         Datam = theta';
%         figure, imshow(tmpim), title(var(Datam));
%         continue;
        if var(sum(Datam))>100
            [cluster_idx, cluster_center] = kmeans(Datam,3,'Start',kMeansClusterCenters);
            
            
            
            pixel_labels = reshape(cluster_idx,mrows,mcols);
            %         figure, imshow(pixel_labels,[]), title('image labeled by cluster index');
            
            for i = 1: mrows
                for j = 1: mcols
                    if (pixel_labels(i, j) == lm_id)
                        pixel_labels(i, j)= 3;
                    else
                        if (pixel_labels(i, j) == nc_id)
                            pixel_labels(i, j)= 1;
                        else
                            if (pixel_labels(i, j) == st_id)
                                pixel_labels(i, j)= 2;
                            end
                        end
                    end
                end
            end
            
        else
            pixel_labels=ones(mrows,mcols)*lm_id;
        end
        
        firstLine = zeros(1,mcols)-1;
        firstLine(1,1)= mrows;
        firstLine(1,2)= mcols;
        KMap = [firstLine ; pixel_labels];
        KmeansImFileName = [fullFileName(indx+1:end-4) '_ImX' num2str(testIndX) '_ImY' num2str(testIndY) '_stX' num2str(ii) '_stY' num2str(jj) '_k3'];
        save([pwd '/' KmeansImFileName],'KMap');
        dlmwrite([pwd '/' KmeansImFileName],KMap);
        testIndY=testIndY+1;
        clear aa b c d;
        %         pack
    end
    testIndY=1;
    testIndX = testIndX+1;
    
end
