clear
% this script is used to calculate the distance between two close membranes: membrane1 and membrane2 (e.g. double membrane of some organellels)
% the output file 'memb_dist_display.bild' shows double-membrane distances distribution can be read in Chimera software.
%
% Notes: this script relies on TOM software toolbox, which contains essential function (e.i. tom_mrcread)
% the package of TOM software toolbox can be fonud:
% Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
% Journal of Structural Biology, 149 (2005), 227-234.

%% color code
rgb1=[216 40 69] ./[255,255,255]; % represent close
rgb2=[0, 189, 242]./[255,255,255]; % represent far

%% range of distance (pixels) between memb1 and memb2
dist_range=[4,16]; % the real distance values should be within this range

%% read tomo
mb1=tom_mrcread('053_memb1_26.64Apx.mrc'); % mrc file with membrane1 segmentation, pixel value of membrane = 1
mb1=mb1.Value;
mb2=tom_mrcread('053_memb2_26.64Apx.mrc'); % mrc file with membrane2 segmentation, pixel value of membrane = 1
mb2=mb2.Value;

%% output parameters
sphere_r=0.5; % the radius of sphere composing the output membrane


fprintf('counting the num of slices with two membranes\n');
zlist=[];
for z=1:size(mb1,3)
    slice1=mb1(:,:,z); % get slice of tomo
    slice2=mb2(:,:,z);
    if ~isempty(find(slice1==1,1)) && ~isempty(find(slice2==1,1))
    zlist=[zlist,z]; % slices contain two memb
    end
end

AEP1=[]; %all edge points
AEP2=[];
EIV1=zeros(size(mb1)); %edge in volume
EIV2=zeros(size(mb2));
fprintf('edge of membrane1 are obtained in every slice, please wait:\n');
for z=zlist
    slice1=mb1(:,:,z); %get slice of tomo
    slice2=mb2(:,:,z);
    fprintf('doing slice%d\n',z);    
    [edgePxy1,edgePxy2]=memb_to_points_in_2Dslice(slice1,slice2,dist_range); %get edge coords
    ePxyz1=[edgePxy1,ones(size(edgePxy1,1),1)*z];
    ePxyz2=[edgePxy2,ones(size(edgePxy2,1),1)*z];
    AEP1=[AEP1;ePxyz1];
    AEP2=[AEP2;ePxyz2];
    for i = 1:size(ePxyz1,1)
        EIV1(ePxyz1(i,1),ePxyz1(i,2),ePxyz1(i,3))=1;
    end
end
fprintf('generating edge of membrane1: edge_vol.mrc \n');
tom_mrcwrite(EIV1,'name','edge_vol.mrc');% output mrc file contains edge of membrane1

D1 = bwdist(mb1); %dist function
D2 = bwdist(mb2);
fprintf('calculating the closest distance from the edge of membrane1 to membrane2\n');
for i = 1:size(AEP1,1)
    x1=AEP1(i,1);
    y1=AEP1(i,2);
    z1=AEP1(i,3);
    dist=D2(x1,y1,z1); %dist between edge
    AEP1(i,4)=dist;
end
save('all_edge_points_dist.mat','AEP1'); % output file contains coordinates and distances
distR=[min(AEP1(:,4)),max(AEP1(:,4))]; %distance range

% delPL=round(1+rand(1,size(AEP1,1)/2)*(size(AEP1,1)-1)); 
% AEP1(delPL,:)=[];

fprintf('generating output file: memb_dist_display.bild \n');
fid=fopen('memb_dist_display.bild','a'); % output file shows double-membrane distances distribution can be read in Chimera software.
for i = 1:size(AEP1,1)
    fprintf(fid,'.color %f %f %f\n',rgb1-(AEP1(i,4)-distR(1))/(distR(2)-distR(1))*(rgb1-rgb2));
    fprintf(fid,'.sphere %d %d %d %d\n',AEP1(i,1)/1,AEP1(i,2)/1,AEP1(i,3)/1,sphere_r);
end
fclose(fid);
fprintf('done! \n');


function [edgePxy1,edgePxy2]=memb_to_points_in_2Dslice(slice1,slice2,dist_range)
% to fill the inter-space between memb1 and memb2
% input: slice1 contains memb1; slice2 contains memb2

%% caculate dist between two memb
[x1,y1] = ind2sub(size(slice1),find(slice1==1));
coords1=[x1,y1];
[x2,y2] = ind2sub(size(slice2),find(slice2==1));
coords2=[x2,y2];

dist=pdist2(coords1,coords2);
k3=find(dist<dist_range(2) & dist>dist_range(1)); % get point pairs (dist smaller than range) in memb1 and memb2 respectively
[p1,p2] = ind2sub(size(dist),k3);
pairs=[p1,p2];

edgePxy1=zeros(size(pairs,1),2);
edgePxy2=zeros(size(pairs,1),2);
for i = 1:size(pairs,1)
    p_id1=pairs(i,1);
    p_id2=pairs(i,2);
    p1xy=coords1(p_id1,:); % (x,y) of p1
    p2xy=coords2(p_id2,:); % (x,y) of p2
    a=(p2xy(2)-p1xy(2))/(p2xy(1)-p1xy(1)); % y=ax+b
    if a ~= Inf && a~=0 && a ~= -Inf
        b=p2xy(2)-a * p2xy(1);
        if a>=1 || a<=-1 %斜率大
            if p1xy(2)<p2xy(2)
                yrange= p1xy(2):1:p2xy(2);
            else
                yrange= p1xy(2):-1:p2xy(2);
            end
            thr=[((yrange-b*ones(size(yrange))) ./ (a*ones(size(yrange))))',yrange'];
            thr=round(thr);
        elseif (a<1 && a>0) || (a>-1 && a<0) %斜率小
            if p1xy(1)<p2xy(1)
                xrange= p1xy(1):1:p2xy(1);
            else
                xrange= p1xy(1):-1:p2xy(1);
            end
            thr=[xrange',(a*ones(size(xrange)).*xrange + b*ones(size(xrange)))'];
            thr=round(thr); % pix coords got through by line
        end        
    elseif a == Inf || a== -Inf %线垂直
        if p1xy(2)<p2xy(2)
            yrange= p1xy(2):1:p2xy(2);
        else
            yrange= p1xy(2):-1:p2xy(2);
        end
        thr=[p1xy(1,1)*ones(size(yrange))' , yrange'];
    elseif a == 0 %线水平
        if p1xy(1)<p2xy(1)
            xrange= p1xy(1):1:p2xy(1);
        else
            xrange= p1xy(1):-1:p2xy(1);
        end
        thr=[xrange',p1xy(1,2)*ones(size(xrange))'];
    end
    
    lineV1=zeros(1,size(thr,1));%记直线上的value
    for m = 1:size(thr,1)
        lineV1(m)=slice1(thr(m,1),thr(m,2)) ; %for memb1, line value (ex.1111110000)
    end    
    for n = 1:size(lineV1,2)-1
        if lineV1(1)==1 && lineV1(n)==1 && lineV1(n+1)==0 
            if isempty(find(lineV1((n+1):size(lineV1,2))==1,1)) %排除11100001100，11100001等情况
                edgeP1=[thr(n,1),thr(n,2)];
                break
            else
                edgeP1=[];
                break
            end
        end
        edgeP1=[];
    end
    
    lineV2=zeros(1,size(thr,1)); %记直线上的value
    for m = 1:size(thr,1)
        lineV2(m)=slice2(thr(m,1),thr(m,2)) ; %for memb2, line value (ex.000001111)
    end    
    for n = 1:size(lineV2,2)-1
        if lineV2(1)==0 && lineV2(n)==0 && lineV2(n+1)==1 
            if isempty(find(lineV2((n+1):size(lineV2,2))==0,1)) %排除000011010，0000111等情况
                edgeP2=[thr(n+1,1),thr(n+1,2)];
                break
            else
                edgeP2=[];
                break
            end
        end
        edgeP2=[];
    end  
    
    if ~isempty(edgeP1) && ~isempty(edgeP2)
        edgePxy1(i,:)=edgeP1;
        edgePxy2(i,:)=edgeP2;
    end     
end

edgePxy1=unique(edgePxy1,'row');
edgePxy2=unique(edgePxy2,'row');
edgePxy1(find(edgePxy1==0,1),:)=[]; %除0
edgePxy2(find(edgePxy2==0,1),:)=[];

end
