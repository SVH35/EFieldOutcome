%% E-Field Outcome Measure Extraction Script
% by Sybren Van Hoornweder (sybren.vanhoornweder@uhasselt.be),
% based on code accompagnying SimNIBS

% this code accompagnies the following publication: UNDER REVIEW

%%% requirements
% requires SimNIBS to be installed and added to the MATLAB environment

%%% input
% pathname = path that contains .msh simulation file
% mni_coord = MNI coordinate of grey matter ROI that should be analyzed
% roi_sphere = size of sphere used for ROI analyses, typically 10mm is used
% toplot = can be either 1 (show plot) or 0 (don't show plot)

%%% output
% structure that contains the following fields; EField in ROI and P99,
% analyzed ROI and P99 grey matter volumes, amount of overlap between
% analyzed volumes and EField ratio. 

%%% example
%outcome = EFieldOutcome('D:\Lab\Sybren\sim_030122', [-52.2, -16.4, 57.8], 10, 1);

%% Function
function outcome = EFieldOutcome(pathname, mni_coord, roi_sphere, toplot)

%% get simulation results, load them them
cd(pathname);
mydir = pwd;

if isunix == 1
    idcs = strfind(mydir,'/');
else     idcs = strfind(mydir,'\');
end

cd(mydir(1:idcs(end)-1))
sub_m2m = dir('m2m*');

cd(pathname);
simulation = dir('*.msh');
sub_efieldmesh = simulation.name;
m = mesh_load_gmsh4(sub_efieldmesh);

%% get grey matter mask
cd(mydir(1:idcs(end)-1))
gray_matter = mesh_extract_regions(m, 'region_idx', 2);     % get grey matter mask
subject_coord = mni2subject_coords(mni_coord, sub_m2m.name);
elm_centers = mesh_get_tetrahedron_centers(gray_matter);    % get element centers
elm_vols = mesh_get_tetrahedron_sizes(gray_matter);         % get volumes for averaging

%% GET ROI EFIELD
roi = sqrt(sum(bsxfun(@minus, elm_centers, subject_coord).^2, 2)) < roi_sphere;  % roi mask

% Get field and calculate mean
field_name = 'normE';
field_idx = get_field_idx(gray_matter, field_name, 'elements');
field = gray_matter.element_data{field_idx}.tetdata;
outcome.roi_efield = sum(field(roi) .* elm_vols(roi))/sum(elm_vols(roi));

%% GET P99
s.field_idx = 2;
s.region_idx = 2;
s.datatype = 'tet';
m=mesh_extract_regions(m,'elemtype','tet','region_idx',s.region_idx);
[data, name, scaleLimits, elemsizes, elempos] = get_data_and_scaleLimits(m,s.field_idx,s.datatype,[]);

idx=~isnan(data);
ElemsizesUnsorted = elemsizes;
data=data(idx);
elemsizes=elemsizes(idx);
elempos=elempos(idx,:);

% sort data, get cdf
UnsortedData = data;
[data,idx] = sort(data);
idx_raw = idx;
elemsizes=elemsizes(idx);
elemsizes=cumsum(elemsizes);
elemsizesNormed=elemsizes/elemsizes(end);
elempos=elempos(idx,:);

s.field_name=name;
s.max=max(data);
s.XYZ_max=elempos(data==s.max,:);
idx=find(elemsizesNormed>99/100,1,'first');
if isempty(idx); idx=length(data); end
outcome.P99_efield=data(idx);

% mean and SD of positions (equations weighted for element size)
meanVal=sum(bsxfun(@times,elempos(idx:end,:),elemsizes(idx:end)),1)./...
    repmat(sum(elemsizes(idx:end)),1,3);

N_nonzero=sum(elemsizes(idx:end,:)>0);
scaleFac=(N_nonzero-1)/N_nonzero*sum(elemsizes(idx:end,:));
for j=1:3
    s.XYZstd_perc(j)=sqrt( sum(elemsizes(idx:end).*( elempos(idx:end,j)-meanVal(j) ).^2,1)/scaleFac );
end
s.XYZ_perc=meanVal;

%% GET ROI and P99 mask
idx=find(elemsizesNormed>99/100,1,'first');
idxlist = find(elemsizesNormed>99/100,length(elemsizesNormed));
p99mask = zeros(length(data),1);
s.p99mask(idxlist) = 1;
temp = sortrows([s.p99mask', idx_raw],2);
s.p99mask(:) = temp(:,1)';

for mask_i = 1:length(field)
    if roi(mask_i)==s.p99mask(mask_i)
        if roi(mask_i)==1
            finalmask(mask_i,1) = 1;
        else finalmask(mask_i,1) = 0; end
    else finalmask(mask_i,1) = 0; end
end

s.finalvolume = ElemsizesUnsorted .* finalmask;
outcome.volume_overlap = sum(s.finalvolume);

if isempty(idx); idx=length(data); end
peakvalue=data(idx);

idx=find(data>=99/100*peakvalue,1,'first');
if idx == 1
    outcome.volume_P99 = elemsizes(end);
else  outcome.volume_P99 =elemsizes(end)-elemsizes(idx-1);
end

outcome.volume_ROI = sum(elm_vols(roi));
outcome.EFieldRatio = outcome.roi_efield / outcome.P99_efield;

%% plot figure
if toplot == 1
    patch('Faces',gray_matter.tetrahedra,'Vertices',gray_matter.nodes,'FaceVertexCData',data,...
        'FaceColor','#808080','EdgeColor','#a9a9a9');
    hold on
    patch('Faces',gray_matter.tetrahedra(roi,1:4),'Vertices',gray_matter.nodes,'FaceVertexCData',data,...
        'FaceColor','#177BA2','EdgeColor','#32b9ed');
    hold on
    patch('Faces',gray_matter.tetrahedra(logical(s.p99mask),1:4),'Vertices',gray_matter.nodes,'FaceVertexCData',data,...
        'FaceColor','#C32127','EdgeColor','#fc5157');
    hold on
    patch('Faces',gray_matter.tetrahedra(logical(s.finalvolume),1:4),'Vertices',gray_matter.nodes,'FaceVertexCData',data,...
        'FaceColor','#B498B5','EdgeColor','#8d598f');
else end
end
