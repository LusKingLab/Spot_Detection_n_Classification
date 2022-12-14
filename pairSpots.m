function [pairs_data, conflict_data] = pairSpots(spotData_1, spotData_2, Rmax)
% UNDER CONSTRUCTION
% pair spots from different spotLists based on proximity.
% Runs through all spots in the spotData_1 and
% Looks for the spots ins postData_2 within Rmax: 
% if no spots - leaves unpaired,
% if one - pairs,
% if more than one  - chooses the closest but makes a note that there was an ambiguity 
%
% INPUT:
%  spotData - cell array where each entry is an individual spot data organized as structure with
%  following fields for i-th spot:
%  spotData{ii} 
%   struct with fields:
%      spotPosition: [x, y, z] of the spot position
%      intensityRatio: 1.2747 ratio calculated at spot identification
%      spotIntensity: 2555326 intenisty of the peak area calculated at the spot identification
%
% OUTPUT:
%   pairs_data - table, where each row 'kk' contains info on paired spots as:
%                            pairs_data(kk,:) = [spot_1 ID (i.e.# in spotData_1), spot_2 ID, distance between spots,...
%                                                                     spot_1 Intensity, spot_2 intensity]
%   conflict_data - table containig all spots that have than one neighbor, where each row 'mm' is 
%                                   (mm, :) = [spot_ID,n_neighbours];
%
% Ivan Surovtsev, 2021.12.03

n_spots_1= length(spotData_1);
n_spots_2= length(spotData_2);
Rmax2=Rmax^2;

pairs_data=[];
conflict_data=[];

% spots_ID=1:n_spots_1;

% make position table and vectors
pos_all_1=zeros(n_spots_1,5);
SI_all_1=zeros(n_spots_1,1);
for ii=1:n_spots_1
    pos_all_1(ii,1:4)=[ii,spotData_1{ii}.spotPosition];
    SI_all_1(ii)=spotData_1{ii}.spotIntensity;
end
all_X1=pos_all_1(:,2);
 all_Y1=pos_all_1(:,3);
 all_Z1=pos_all_1(:,4);

pos_all_2=zeros(n_spots_2,5);
SI_all_2=zeros(n_spots_2,1);
for ii=1:n_spots_2
    pos_all_2(ii,1:4)=[ii,spotData_2{ii}.spotPosition];
    SI_all_2(ii)=spotData_2{ii}.spotIntensity;
end
all_X2=pos_all_2(:,2);
 all_Y2=pos_all_2(:,3);
 all_Z2=pos_all_2(:,4);

% start ith the frist spot find the closest negbour in the other spotList 
% update whatever left list and continue...
% n_pairs=0;
% spots_left=spots_ID;

XX1=repmat(all_X1,1,n_spots_2);
 YY1=repmat(all_Y1,1,n_spots_2);
 ZZ1=repmat(all_Z1,1,n_spots_2);
 
XX2=repmat(all_X2',n_spots_1,1);
 YY2=repmat(all_Y2',n_spots_1,1);
 ZZ2=repmat(all_Z2',n_spots_1,1); 

RR2=(XX1-XX2).^2 + (YY1-YY2).^2+ (ZZ1-ZZ2).^2;

ind=RR2<Rmax2;

kk=0;
mm=0;
for ii=1:n_spots_1
    n_neighbours=sum(ind(ii,:));
    if n_neighbours==1
        kk=kk+1;
        [dist,jj]=min(RR2(ii,:));
        pairs_data(kk,:)=[ii,jj,dist,SI_all_1(ii), SI_all_2(jj)];        
    elseif n_neighbours >1
        %currently not solved for potential conflicts
         kk=kk+1;
        mm=mm+1;
        [dist,jj]=min(RR2(ii,:));
        pairs_data(kk,:)=[ii,jj,dist,SI_all_1(ii), SI_all_2(jj)];
        conflict_data(mm,:)=[ii,n_neighbours];
    end
    
end



disp(['paired ',num2str(kk),' spots'])
disp([num2str(mm), ' spots had more than 1 neighbour '])
disp([num2str(n_spots_1-kk), ' and ',num2str(n_spots_2-kk) ' in 1st and 2nd spotLists,respectively, did not have neighbours '])



end