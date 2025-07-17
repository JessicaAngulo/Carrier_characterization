%% Image Acquisition and File Organization
% This code imports all images. They must be in TIFF format (up to 16-bit)
% and the name of the file should not contain points or weird simbols (/,\,*, etc). 
% Ensure that the channel to be quantified is split beforehand.

[filename, path] = uigetfile('*.tif','Select TIFF images','MultiSelect','on');
path = string(path) + string(filename);
n_img = length(filename);

%% Cell Segmentation
% Segmentation of each cell in every image is done manually.
% Drawn ROIs are stored along with the computed area.

n_sets = length(path);
cell_rois = cell(n_sets,1);

for i = 1:n_sets
    info = imfinfo(path{i});
    image_w = info.Width;
    image_h = info.Height;
    BitDepth = info.BitDepth;

    if BitDepth == 8
        image_i = zeros(image_h, image_w, 'uint8');
    elseif BitDepth == 16
        image_i = zeros(image_h, image_w, 'uint16');
    else
        image_i = zeros(image_h, image_w, 'double');
    end

    tstack = Tiff(path{i}, 'r');
    image_i = tstack.read();
    image_i = imadjust(image_i);
    image_i = uint8(image_i / 256);

    imshow(image_i);
    for j = 1:20
        try
            roi = drawpolygon('FaceAlpha',0,'Deletable',0,'InteractionsAllowed','none',...
                'LineWidth',1,'Color','r','DrawingArea',[0,0,image_w,image_h]);
            if length(roi.Position) > 2
                cell_rois{i}{1,j} = roi.Position;
                cell_rois{i}{2,j} = polyarea(roi.Position(:,1), roi.Position(:,2));
            end
        catch
            break
        end
    end
    log = cellfun(@isempty,cell_rois{i}(1,:));
    cell_rois{i}(:,log) = [];
end

%% Threshold Selection
% Choose a representative image and define a threshold for vesicle selection.
% Same threshold will be used for all fields.

i = 1; % Change to image index of interest
info = imfinfo(path{i});
image_w = info.Width;
image_h = info.Height;
BitDepth = info.BitDepth;

if BitDepth == 8
    image = zeros(image_h, image_w, 'uint8');
elseif BitDepth == 16
    image = zeros(image_h, image_w, 'uint16');
else
    image = zeros(image_h, image_w, 'double');
end

tstack = Tiff(path{i}, 'r');
image = tstack.read();

image_j = imgaussfilt(image,1);
se = strel("sphere",1);
image_j = imerode(image_j, se);
image_j = uint8(image_j / 256);

vertices = cell(1, width(cell_rois{i}));
for j = 1:width(cell_rois{i})
    vertices{j} = cell_rois{i}{1,j};
end
binary_roi = poly2label(vertices, 1:width(vertices), [image_h, image_w]);
image_j(~binary_roi) = 0;

th = [mean(image_j,'all') + std(double(image_j), 0, 'all'), 256];
binary = image_j >= th(1) & image_j <= th(2);

% Display threshold result
redChannel = image_j;
greenChannel = image_j;
blueChannel = image_j;
redChannel(binary) = 256;
greenChannel(binary) = 0;
blueChannel(binary) = 0;
rgbImage = cat(3, redChannel, greenChannel, blueChannel);
imshow(rgbImage);

%% Vesicle Segmentation
% Segments structures within a defined size range and extracts features:
% Area in μm², Centroid and Median brightness intensity

min_s = 5;
max_s = 500;
pixel_size = 0.06;
boundaries = cell(1,n_sets);
centroids = cell(1,n_sets);
areas = cell(1,n_sets);
brightness_intensity = cell(1,n_sets);
median_param = NaN(n_sets, 4);

for i = 1:n_sets
    info = imfinfo(path{i});
    image_w = info.Width;
    image_h = info.Height;
    BitDepth = info.BitDepth;

    if BitDepth == 8
        image = zeros(image_h, image_w, 'uint8');
    elseif BitDepth == 16
        image = zeros(image_h, image_w, 'uint16');
    else
        image = zeros(image_h, image_w, 'double');
    end

    tstack = Tiff(path{i}, 'r');
    image = tstack.read();

    image_j = imgaussfilt(image, 1);
    se = strel("sphere", 1);
    image_j = imerode(image_j, se);
    image_j = uint8(image_j / 256);

    vertices = cell(1, width(cell_rois{i}));
    for r = 1:width(cell_rois{i})
        vertices{r} = cell_rois{i}{1,r};
    end
    binary_roi = poly2label(vertices, 1:width(vertices), [image_h, image_w]);
    image_j(~binary_roi) = 0;

    th = [mean(image_j, 'all') + std(double(image_j), 0, 'all'), 256];
    binary = image_j >= th(1) & image_j <= th(2);
    image(~binary) = 0;

    boundary = bwboundaries(binary, 'noholes');
    boundaries_j = {};
    areas_j = [];
    centroids_j = [];
    brightness_j = [];

    for k = 1:length(boundary)
        area = polyarea(boundary{k}(:,1), boundary{k}(:,2));
        if area > min_s && area < max_s
            polyin = polyshape(boundary{k});
            [x, y] = centroid(polyin);
            center = [x, y];
            mask = poly2mask(boundary{k}(:,2), boundary{k}(:,1), image_w, image_h);
            b_intensity = median(double(image(mask)));

            boundaries_j = [boundaries_j; boundary{k}];
            areas_j = [areas_j; area * pixel_size^2];
            centroids_j = [centroids_j; center];
            brightness_j = [brightness_j; b_intensity];
        end
    end

    boundaries{i} = boundaries_j;
    areas{i} = areas_j;
    centroids{i} = centroids_j;
    brightness_intensity{i} = brightness_j;
    median_param(i,1:3) = [length(areas_j), median(areas_j), median(brightness_j)];
end

% Record cell area
for i = 1:n_sets
    if ~isempty(cell_rois{i})
        median_param(i,4) = cell_rois{i}{2,1} * pixel_size^2;
    end
end

%% Visualization
% Overlays segmented vesicles on the image using label overlays.

i = 1; % Select image for visualization
info = imfinfo(path{i});
image_w = info.Width;
image_h = info.Height;
BitDepth = info.BitDepth;

if BitDepth == 8
    image = zeros(image_h, image_w, 'uint8');
elseif BitDepth == 16
    image = zeros(image_h, image_w, 'uint16');
else
    image = zeros(image_h, image_w, 'double');
end

tstack = Tiff(path{i}, 'r');
image = tstack.read();
image = uint8(image / 256);

vertices = cell(1, height(boundaries{i}));
for j = 1:height(boundaries{i})
    vertices{j} = boundaries{i}{j};
end
L = poly2label(vertices, 1:width(vertices), [image_w, image_h]);
I2 = labeloverlay(imadjust(image), L, 'Colormap', 'hsv');
imshow(I2);

