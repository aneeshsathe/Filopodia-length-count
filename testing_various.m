clear
clc
%%
% addpath('Z:\IT Core Work\Naila_filopodia\modules\')

add_path1=fullfile(fileparts(mfilename('fullpath')),'modules', 'image_tools', 'bfmatlab', filesep);
add_path2=fullfile(fileparts(mfilename('fullpath')),'modules', 'dlmcell', filesep);
add_path3=fullfile(fileparts(mfilename('fullpath')), 'modules', filesep);
addpath(add_path1)
addpath(add_path2)
addpath(add_path3)
%%

file_path='Z:\IT Core Work\Naila_filopodia\140826_HeLaJW_GFPmyo10_CheLifeA  5 0.03tps with and wo ATR drug\140826_HeLaJW_GFPmyo10_CheLA  5 0.03tps.tif';

img_write_path='Z:\IT Core Work\Naila_filopodia\control_';
%%
[img1,img2, Rimg1, Rimg2]=read_im(file_path);
% [img1,img2 Rimg1, Rimg2]=read_im();
%% Ratio

rat_img=Rimg1./Rimg2;

bwrat_img=im2bw(rat_img,graythresh(rat_img));
mImg2=bwrat_img.*img2;
mImg1=bwrat_img.*img1;
% imtool([mImg1 img1 mImg2 img2 mImg1+mImg2],[])

MCimg=mat2gray(mImg1+mImg2); %add channels
[cell_body, unerod_cell_body]=keep_max_area_obj(an_gray(MCimg));
cell_perim=bwperim(cell_body);

% imshow(cell_perim)
%% texture based

bw_tex=texture_analy(MCimg );

% imshow(cell_body)
some_stuff=(bw_tex.*MCimg)+MCimg-cell_body;
% imtool(some_stuff+cell_perim)

se = strel('diamond', 4);
aft_open=imdilate(cell_perim,se);
remove_perim=an_gray(some_stuff).*imcomplement(aft_open);

% imshow(aft_open)
% imshow([bw_tex some_stuff remove_perim])
% +cell_perim MCimg
% area thresh
%%
out_BW = area_thresh_obj( remove_perim, 4 );

% imshow((out_BW.*imcomplement(cell_body))+cell_perim)
comb_img=img1+img2;
%%
% imshow([img1, img2, img1+img2 ; cell_body, out_BW.*imcomplement(cell_body), out_BW.*imcomplement(cell_body)+cell_perim  ] )
%%

filoBW_img=out_BW.*imcomplement(cell_body);
%%
% imtool([out_BW, filoBW_img.*MCimg])
filo_prop=labelmatrix(bwconncomp(filoBW_img));
% imtool(filo_prop)


% imtool(label2rgb(filo_prop,'spring', 'c', 'shuffle'))

% figure,
% imshow(MCimg)

[xls_label{1,1:3}]=deal('No.','ChessDist','EucDist');


for count=1:max(max(filo_prop))
    %%
    bw2=zeros(size(filo_prop));
    % bw1(filo_prop==2)=comb_img(filo_prop==2);
    bw1=(filo_prop==count);
    
    
    f=regionprops(bw1,'Extrema');
    g=regionprops(cell_perim,'Centroid','Area');
    cell_cent=floor(g.Centroid);
    b=ceil(f.Extrema);
    bw2(b(:,2),b(:,1))=5;
    bw2(cell_cent(2),cell_cent(1))=6;
    %  imtool(bw2+bw1+cell_perim)
    %%
    filo_skel = bwmorph(bw1, 'skel', Inf);
    %     imshow(filo_skel)
    
    %%
    % imtool(bwlabel(bw1+cell_perim))
    D1 = bwdist(bw1, 'cityblock');
    D2 = bwdist(cell_perim, 'cityblock');
    D_sum = D1 + D2;
    path_pixels = imregionalmin(D_sum);
    % imshow(path_pixels)
    
    %%
    % [ d,p,a,b ] = an_distofbw( filo_skel,cell_cent,'max' )
    filo_skel_far_pt  = an_nearfar_pt( filo_skel,cell_cent,'max' );
    path_near_pt=an_nearfar_pt( path_pixels,cell_cent,'min' );
    
    % [a,b]=find(filo_skel==1);
    % [d,p] = max(sqrt((a-cell_cent(2)).^2 + (b-cell_cent(1)).^2));
    % d
    %  d = sqrt(d)
    %%
    coll_img(:,:,1)=mat2gray(filo_skel+path_pixels);
    coll_img(:,:,2)=mat2gray(cell_perim+path_pixels);
    coll_img(:,:,3)=mat2gray(bw1+path_pixels);
    %     imagesc(coll_img)
    %     hold on
    %     plot(filo_skel_far_pt(1), filo_skel_far_pt(2), 'g*', 'MarkerSize', 5)
    %     plot(path_near_pt(1), path_near_pt(2), 'r*', 'MarkerSize', 5)
    %     hold off
    %%
    
    all_img=im2bw(filo_skel+cell_perim+path_pixels);
    % %     imshow(all_img)
    %     hold on
    %     plot(filo_skel_far_pt(1), filo_skel_far_pt(2), 'g*', 'MarkerSize', 5)
    %     plot(path_near_pt(1), path_near_pt(2), 'r*', 'MarkerSize', 5)
    %     hold off
    
    %%
    %     D1 = bwdistgeodesic(all_img, filo_skel_far_pt(1), filo_skel_far_pt(2), 'quasi-euclidean');
    %     D2 = bwdistgeodesic(all_img, path_near_pt(1), path_near_pt(2), 'quasi-euclidean');
    D1 = bwdistgeodesic(all_img, filo_skel_far_pt(1), filo_skel_far_pt(2),'chessboard');
    D2=bwdistgeodesic(all_img, path_near_pt(1), path_near_pt(2),'chessboard');
    
    
    
    D = D1 + D2;
    D = round(D * 8) / 8;
    
    D(isnan(D)) = inf;
    skeleton_path = imregionalmin(D);
    
    path_length = D(skeleton_path);
    path_length = path_length(1);
    chess_dist=path_length;
    %%
    % P = imoverlay(all_img+comb_img, imdilate(skeleton_path, ones(1,1)), [1 0 0]);
    %         figure, imshow(P, 'InitialMagnification', 200)
    P1 = imoverlay(imoverlay([comb_img,comb_img], skeleton_path, [1 0 0]),cell_perim,[0,1,0]);
    
    %     out_P=insertText(P,[path_near_pt(1)+5, path_near_pt(2)+5],['Dist:', num2str(path_length), 'px']);
    out_P1=insertText(P1,[1,1],['No. ',num2str(count), ' | Dist:', num2str(path_length), ' px. | Method: Chessboard']);
    
    pos   = [filo_skel_far_pt(1), filo_skel_far_pt(2);path_near_pt(1), path_near_pt(2);...
        filo_skel_far_pt(1)+512, filo_skel_far_pt(2);path_near_pt(1)+512, path_near_pt(2)];
    color = {'yellow', 'yellow','yellow', 'yellow'};
    out_P1 = insertMarker(out_P1, pos, 's', 'color', color, 'size', 3);
    
    %%
    D1 = bwdistgeodesic(all_img, filo_skel_far_pt(1), filo_skel_far_pt(2), 'quasi-euclidean');
    D2 = bwdistgeodesic(all_img, path_near_pt(1), path_near_pt(2), 'quasi-euclidean');
    
    
    
    
    D = D1 + D2;
    D = round(D * 8) / 8;
    
    D(isnan(D)) = inf;
    skeleton_path = imregionalmin(D);
    
    path_length = D(skeleton_path);
    path_length = path_length(1);
    eu_dist=path_length;
    %%
    % P = imoverlay(all_img+comb_img, imdilate(skeleton_path, ones(1,1)), [1 0 0]);
    %         figure, imshow(P, 'InitialMagnification', 200)
    P2 = imoverlay(imoverlay([comb_img,comb_img], skeleton_path, [1 0 0]),cell_perim,[0,1,0]);
    
    %     out_P=insertText(P,[path_near_pt(1)+5, path_near_pt(2)+5],['Dist:', num2str(path_length), 'px']);
    out_P2=insertText(P2,[1,1],['No. ',num2str(count), ' | Dist:', num2str(path_length), ' px. | Method: quasi-euclidean']);
    
    pos   = [filo_skel_far_pt(1), filo_skel_far_pt(2);path_near_pt(1), path_near_pt(2);...
        filo_skel_far_pt(1)+512, filo_skel_far_pt(2);path_near_pt(1)+512, path_near_pt(2)];
    color = {'yellow', 'yellow','yellow', 'yellow'};
    out_P2 = insertMarker(out_P2, pos, 's', 'color', color, 'size', 3);
    %     imshow(out_P)
    %% Collect Data
    [pro_xls_label{1,1:3}]=deal(['No ',num2str(count)],chess_dist,eu_dist);
    xls_label=vertcat(xls_label,pro_xls_label);
%     [pro_xls_label{1,1:3,count}]=deal(['No ',num2str(count)],chess_dist,eu_dist);
    
    
    
%     xls_label=vertcat(xls_label,reshape(permute(pro_xls_label, [1 3 2]),[],3));
    %% Write Image
    out_P=[out_P1;out_P2];
    if count==1
        imwrite(out_P,[img_write_path,'out_img.tif'])
    else
        imwrite(out_P,[img_write_path,'out_img.tif'],'WriteMode','append')
    end
    disp(['No.',num2str(count),' of ', num2str(max(max(filo_prop)))])
    %     imshow(out_P);
    %     axis off
    %     axis image
    % %     text(path_near_pt(1)+5, path_near_pt(2)+5,['Dist:', num2str(path_length), 'px'],...
    % %         'BackgroundColor',[.7 .9 .7],...
    % %         'HorizontalAlignment','left','VerticalAlignment','Top',...
    % %         'FontSize',7)
    %     hold on
    %     plot(filo_skel_far_pt(1), filo_skel_far_pt(2), 'y.', 'MarkerSize', 15)
    %     plot(path_near_pt(1), path_near_pt(2), 'y.', 'MarkerSize', 15)
    %     hold off
    %
end

%% write data

[xls_label{1,4:5}]=deal('Perim ',g.Area);
%%
xls_file_name=[img_write_path,'out_results.txt'];
    
    dlmcell(xls_file_name,xls_label,',')
   
    %write excel file if OS is windows
    if ispc
        xls_file_name2=[img_write_path,  'out_results.xls'];
        xlswrite(xls_file_name2,xls_label,'Sheet1' );
        xlswrite(xls_file_name2,xls_label,'Sheet2' );
    end
    
    
    disp('done writing excel file')
%%
% imshow(mat2gray(cell_perim+skeleton_path))
%
% hold on
% plot(filo_skel_far_pt(1), filo_skel_far_pt(2), 'go', 'MarkerSize', 5)
% plot(path_near_pt(1), path_near_pt(2), 'bo', 'MarkerSize', 5)
% hold off





%% Homomorphic filtering

% [ Ihmf_3, Ihmf_2, Ihmf_spatial, Ihmf ] = homomorph_filt(MCimg);

% imtool([ Ihmf_3 Ihmf_2 Ihmf_spatial Ihmf ])

% imtool(some_stuff+keep_max_area_obj(an_gray(img2))-keep_max_area_obj(an_gray( img1)),[])
%  imshow(an_gray(Ihmf_3))


%% finding shortest paths
% http://blogs.mathworks.com/steve/2011/12/13/exploring-shortest-paths-part-5/
