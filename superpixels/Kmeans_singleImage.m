function Kmeans_singleImage(fullFileName,k)

temp_im = imread(fullFileName);

nrows  = size(temp_im,1);
ncols= size(temp_im,2);
%RGB
Data = double([reshape(temp_im(:,:,1),nrows*ncols,1) reshape(temp_im(:,:,2),nrows*ncols,1) reshape(temp_im(:,:,3),nrows*ncols,1)]);

kMeansClusterCenters = initializeClusterVectorsUsingPrincipalComponent(Data,k);


% determine lumen,stroma, nuclei labels RGB
m= kMeansClusterCenters;
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

KmeansImFileName = [fullFileName(indx+1:end-4) '_k3'];

[cluster_idx, cluster_center] = kmeans(Data,k,'Start',kMeansClusterCenters);

pixel_labels = reshape(cluster_idx,nrows,ncols);

for i = 1: nrows
    for j = 1: ncols
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

firstLine = zeros(1,ncols)-1;
firstLine(1,1)= nrows;
firstLine(1,2)= ncols;
KMap = [firstLine ; pixel_labels];
save([pwd '/' KmeansImFileName],'KMap');
dlmwrite([pwd '/' KmeansImFileName],KMap);