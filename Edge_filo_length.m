clear
clc
%http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/#edgelink
file_path='C:\Users\Nyx\Desktop\Naila_filopodia\Filopodia round2\140826_XY point 14\';
%%
img_write_path=file_path;
Rimg1=imread([file_path,'140826 XY point 14 464ATR 1h17min.tif'],'Index',1);%filo
Rimg2=imread([file_path,'140826 XY point 14 464ATR 1h17min.tif'],'Index',2);%actin

img1=mat2gray(Rimg1);%grey filo
img2=mat2gray(Rimg2);%grey actin

%% crop
imshow(img1+img2,'InitialMagnification',200)
h=imfreehand;
BW_crop=createMask(h);
%
% % Mimg1=Mimg1.*BW_crop;
% % imshow(Mimg1)

%%


Mimg1=mat2gray(adapthisteq(Rimg1)+adapthisteq(Rimg2));
% img1=mat2gray(Rimg1+Rimg2);



% imshow(bwmorph(an_gray(Mimg1),'tophat'))

%%
rat_img=(Rimg1.*2)./Rimg2;

bwrat_img=im2bw(rat_img,graythresh(rat_img));
mImg2=bwrat_img.*img2;
mImg1=bwrat_img.*Mimg1;
% imtool([mImg1 img1 mImg2 img2 mImg1+mImg2],[])

MCimg=mat2gray(mImg1+(mImg2.*2)); %add channels
[cell_body, unerod_cell_body]=keep_max_area_obj(an_gray(MCimg));
cell_perim=bwperim(cell_body);
se = strel('diamond', 4);
aft_open=imdilate(cell_perim,se);

% imshow(aft_open)
%%

INITPSF= fspecial('gaussian',5,5);


[J, P1] = deconvblind(imsharpen(Mimg1,'Radius',2,'Amount',5),INITPSF,5);
[J2, P2] = deconvblind(J,P1,5);

J2=bwareaopen(an_gray(wiener2(J2)),10);
% luc1 = deconvlucy(img1,INITPSF,5);
J3=bwareaopen(an_gray(J2+Mimg1*2),30);
% imtool([J3,Mimg1])
% imshow(an_gray(J2),[])

%%
BW3 = bwmorph((bwmorph(J3-cell_body,'thin',2)+cell_perim).*double(BW_crop),'skel',Inf);
BW4 = bwmorph(bwmorph(bwmorph(bwmorph(BW3,'clean'),'spur'),'clean'),'spur');
% imshow(BW4)
% figure,
% imshow(imoverlay(BW3,BW4,[1 0 0]))
%%

BWbranchpt=bwmorph(BW4,'branchpoints');
BWbranchendpt=bwmorph(BW4,'endpoints');
%%
BW_sans_branchpt=BW4-BWbranchpt;
BW_sans_branchpt_area_filt = bwmorph(bwmorph(bwareaopen(BW_sans_branchpt, 5,8),'clean'),'spur');
% imshow([BW_sans_branchpt_area_filt+BWbranchpt BW_sans_branchpt BW4])
% imshow(BW4-BWbranchendpt)

%%

% disp_img=imoverlay(imoverlay(imoverlay(imoverlay(Mimg1,BW3,[1 0 0]),cell_perim,[0 1 0]),BWbranchpt,[0 0 1]),BWbranchendpt,[1 1 0]);
% imshow(disp_img,'InitialMagnification',200)

%% pixels commom to branch end points and cell_perim

BW_perim_and_branch=BWbranchpt & (cell_perim-BWbranchendpt);

[perim_R,perim_C]=find(BW_perim_and_branch==1);
[branch_R,branch_C]=find(BWbranchendpt==1);

%%

% figure,
% imshow(imoverlay(Mimg1, BW4,[1 0 0]))
% hold on
% for pl_count=1:size(perim_R,1)
%     scatter(perim_C(pl_count),perim_R(pl_count),'*y')
%     text(perim_C(pl_count),perim_R(pl_count),num2str(pl_count))
% end
% for pl_count=1:size(branch_C,1)
%     scatter(branch_C(pl_count),branch_R(pl_count),'*c')
%     text(branch_C(pl_count),branch_R(pl_count),num2str(pl_count))
% end
% hold off
%% Calculate all D1 and D2
% [D1,D2]=deal(NaN(size(BW4,1),size(BW4,2),size(perim_C,1)),NaN(size(BW4,1),size(BW4,2),size(branch_C,1)));
% disp('Calculating Distance Matrices.')
% parfor count_perim =1:size(perim_C,1)
%     pro_D1(:,:,count_perim) = an_D1_D2_find( BW4, perim_C(count_perim), perim_R(count_perim) );
% end
% disp('D1 Done.')
% disp('Calculating D2')
% parfor count=1:size(branch_C,1)
%     pro_D2(:,:,count) = an_D1_D2_find( BW4, branch_C(count), branch_R(count) );
% end
% disp('Dist Matrices Done')

%%
% count=5;
[xls_label{1,1:3}]=deal('No.','ChessDist','EucDist');
tic
for count=1:size(branch_C,1)
    % count_perim=134;
%     in_D2=pro_D2(:,:,count);
%     path_length=zeros(size(perim_C,1),1);
    
        BR_C=branch_C(count);
        BR_R=branch_R(count);
    %     count_perim=size(perim_C,1);
%     disp(['count' num2str(count)])
    parfor count_perim =1:size(perim_C,1)

        
% path_length(count_perim,1)=an_dist_find_filo( pro_D1(:,:,count_perim),in_D2)
        %         pro_path_length= an_dist_find_filo( BW4,perim_C(count_perim), perim_R(count_perim),BR_C, BR_R);
                path_length(count_perim,1)=an_dist_find_filo( BW4,perim_C(count_perim), perim_R(count_perim),BR_C, BR_R)
        
        %         path_length(count_perim)=pro_path_length;
        %         clear pro_path_length
        %         D1 = bwdistgeodesic(BW4, perim_C(count_perim), perim_R(count_perim),'chessboard');
        %         D2=bwdistgeodesic(BW4, branch_C(count), branch_R(count),'chessboard');
        %
        %         %     imshow(cell_perim)
        %         %     hold on
        %         %     scatter(branch_C(1200), branch_R(1200))
        %         D = D1 + D2;
        %         D = round(D * 8) / 8;
        %
        %         D(isnan(D)) = inf;
        %         skeleton_path = imregionalmin(D);
        %
        %         pro_path_length = D(skeleton_path);
        %         path_length(count_perim) = pro_path_length(1);
        %         chess_dist=path_length;
    end
 
    %%
    [path_length, count_perim]=min(path_length);
%     disp(['path length',num2str(path_length)])
    if path_length>5&&path_length<Inf
        D1 = bwdistgeodesic(BW4, perim_C(count_perim), perim_R(count_perim),'chessboard');
        D2=bwdistgeodesic(BW4, branch_C(count), branch_R(count),'chessboard');
        
        %     imshow(cell_perim)
        %     hold on
        %     scatter(branch_C(1200), branch_R(1200))
        D = D1 + D2;
        D = round(D * 8) / 8;
        
        D(isnan(D)) = inf;
        skeleton_path = imregionalmin(D);
        path_length_out = D(skeleton_path);
        path_length_out = path_length_out(1);
        chess_dist=path_length_out;
        % count_perim=min(path_length)
        % P = imoverlay(all_img+comb_img, imdilate(skeleton_path, ones(1,1)), [1 0 0]);
        %         figure, imshow(P, 'InitialMagnification', 200)
        P1 = imoverlay(imoverlay([Mimg1,Mimg1], skeleton_path, [1 0 0]),cell_perim,[0,1,0]);
        
        %     out_P=insertText(P,[path_near_pt(1)+5, path_near_pt(2)+5],['Dist:', num2str(path_length), 'px']);
        out_P1=insertText(P1,[1,1],['No. ',num2str(count), ' | Dist:', num2str(path_length), ' px. | Method: Cityblock']);
        
        pos   = [branch_C(count),  branch_R(count);perim_C(count_perim), perim_R(count_perim);...
            branch_C(count)+512, branch_R(count);perim_C(count_perim)+512, perim_R(count_perim)];
        %     color = {'yellow', 'yellow','yellow', 'yellow'};
        color = {'red', 'red','red', 'red'};
        out_P1 = insertMarker(out_P1, pos, 's', 'color', color, 'size', 3);
        
        %%%%%%%%%%%%%%%%%%%%%
        D1 = bwdistgeodesic(BW4, perim_C(count_perim), perim_R(count_perim),'quasi-euclidean');
        D2 = bwdistgeodesic(BW4, branch_C(count), branch_R(count), 'quasi-euclidean');
        
        D = D1 + D2;
        D = round(D * 8) / 8;
        
        D(isnan(D)) = inf;
        skeleton_path = imregionalmin(D);
        
        path_length_out = D(skeleton_path);
        path_length_out = path_length_out(1);
        eu_dist=path_length_out;
        %%
        % P = imoverlay(all_img+comb_img, imdilate(skeleton_path, ones(1,1)), [1 0 0]);
        %         figure, imshow(P, 'InitialMagnification', 200)
        P2 = imoverlay(imoverlay([Mimg1,Mimg1], skeleton_path, [1 0 0]),cell_perim,[0,1,0]);
        
        %     out_P=insertText(P,[path_near_pt(1)+5, path_near_pt(2)+5],['Dist:', num2str(path_length), 'px']);
        out_P2=insertText(P2,[1,1],['No. ',num2str(count), ' | Dist:', num2str(path_length), ' px. | Method: Cityblock']);
        
        pos   = [branch_C(count),  branch_R(count);perim_C(count_perim), perim_R(count_perim);...
            branch_C(count)+512, branch_R(count);perim_C(count_perim)+512, perim_R(count_perim)];
        %     color = {'yellow', 'yellow','yellow', 'yellow'};
        color = {'red', 'red','red', 'red'};
        out_P2 = insertMarker(out_P2, pos, 's', 'color', color, 'size', 3);
        %     imshow(out_P)
        %% Collect Data
        [pro_xls_label{1,1:3}]=deal(['No ',num2str(count)],chess_dist,eu_dist);
        xls_label=vertcat(xls_label,pro_xls_label);
        
        %% Write Image
        out_P=[out_P1;out_P2];
        if count==1
            imwrite(out_P,[img_write_path,'out_img.tif'])
        else
            imwrite(out_P,[img_write_path,'out_img.tif'],'WriteMode','append')
        end
        disp(['No.',num2str(count),' of ', num2str(size(branch_C,1))])
    end
end
toc
%% write data
g=regionprops(cell_perim,'Area');
[xls_label{1,4:5}]=deal('Perim ',g.Area);
%
xls_file_name=[img_write_path,'out_results.txt'];

dlmcell(xls_file_name,xls_label,',')

%write excel file if OS is windows
if ispc
    xls_file_name2=[img_write_path,  'out_results.xls'];
    xlswrite(xls_file_name2,xls_label,'Sheet1' );
    xlswrite(xls_file_name2,xls_label,'Sheet2' );
end


disp('done writing excel file')
% figure,
% imshow(out_P1)
