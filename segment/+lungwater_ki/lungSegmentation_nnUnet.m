function [mask, lungVolume] = lungSegmentation_nnUnet(image, resolution, tf)

if nargin<3
    tf=1; %default lung segmentation only performed in first timeframe
end

% image_array = permute(abs(image(:,:,:,tf)), [3 4 2 1]); 
image_array = permute(image, [4 3 2 1]);
pxs_spacing = [10 resolution(1) resolution(2)];

% JSON Path with Preprocess Info
path_info_preprocess ='info_preprocess_Dataset011.json';

% Trained Lung Segmentation Model Path
path_onnx = 'GPU_UNet2D_Dataset011_Lung_final_fold1.onnx';

% Preload onnx Network and info_preprocess 
info_preprocess = load_json(path_info_preprocess);
network{1} =importNetworkFromONNX(path_onnx); %PD

disp('Lung segmentation may take up to 3 minutes on CPU.')
mask = eval_case_onnx_numpy_only(image_array, pxs_spacing, network, info_preprocess);
mask = single(permute(mask, [4 3 2 1])); 
lungVolume=sum(mask(:))*resolution(1)*resolution(2)*resolution(3)/10^6; %in ml

%     % Display Mask Details
%     disp(['Mask Ones Count: ', num2str(sum(mask(:) == 1)), '  ', 'Mask Zeros Count: ', num2str(sum(mask(:) == 0))])
%     lungwater_ki.LWDMap_overlay(single(mask_sq),single(mask_sq),data,'Lung mask',[], 0, 1, 0, 1);
%     d(i)=dice(double(squeeze(matfile.setstruct(1).LungWater.LungMask)),double(squeeze(mask_sq)));

end

%--------------------------------------------------------------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------------------------------------------%
% Main Function 
%--------------------------------------------------------------------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------------------------------------------------------------------%
%% Don't change anything from here


%% Main Segmentation Function
function prediction = eval_case_onnx_numpy_test(image_array, pxs_spacing, path_onnx, path_info_preprocess) %loads network and performs segmentation
    info_preprocess = load_json(path_info_preprocess);

    %network = onnxruntime.InferenceSession(path_onnx);
    network = {}; 
    network{end+1} = importONNXNetwork(path_onnx,'OutputLayerType','pixelclassification');

    [data_test, data_properties] = prepare_case_onnx_numpy(image_array, pxs_spacing, info_preprocess);
    data_predicted = predict_case_onnx(data_test, info_preprocess, network);
    prediction = export_prediction_from_softmax_onnx(data_predicted, info_preprocess, data_properties);
end


function prediction = eval_case_onnx_numpy_only(image_array, pxs_spacing, network, info_preprocess) %assumes network has already been loaded
    
    [data_test, data_properties] = prepare_case_onnx_numpy(image_array, pxs_spacing, info_preprocess);
    
    data_predicted = predict_case_onnx(data_test, info_preprocess, network);
    
    prediction = export_prediction_from_softmax_onnx(data_predicted, info_preprocess, data_properties);
    
end


%%

function sigmoid_output = np_sigmoid(x)
    sigmoid_output = 1 ./ (1 + exp(-x));
end


%%
function a = load_json(filepath)
    fid = fopen(filepath, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    a = jsondecode(str);
end

%%
function prediction = predict_case_onnx(data_test, info_preprocess, network) 
    prediction = predict_sliding_window_return_logits(network, data_test, info_preprocess.num_seg_heads, ...
                                                      info_preprocess.patch_size, ...
                                                      'mirror_axes', info_preprocess.inference_allowed_mirroring_axes, ...
                                                      'tile_step_size', info_preprocess.tile_step_size, ...
                                                      'use_gaussian', info_preprocess.use_gaussian, ...
                                                      'precomputed_gaussian', info_preprocess.inference_gaussian);
end

%%
function [data, data_properties] = prepare_case_onnx_numpy(data, pxs_spacing, info_preprocess)
    transpose_f = info_preprocess.transpose_forward;
    IntProp = info_preprocess.foreground_intensity_properties_per_channel;
    target_spacing = info_preprocess.target_spacing;
    normalization_schemes = info_preprocess.normalization_schemes;
    use_mask_for_norm = info_preprocess.use_mask_for_norm;
    
    % Transpose the data
    transpose_f = transpose(transpose_f);
    transpose_vector = [1, arrayfun(@(x) x + 2, transpose_f)];
    data = permute(data, transpose_vector);
    
    data_shape = size(data);
    
    % Extract original spacing
    transpose_f = transpose_f + 1;
    original_spacing = pxs_spacing(transpose_f);
    % Shape before cropping
    shape_before_cropping = data_shape(2:end);
    % Crop to nonzero
    [data, seg, bbox] = crop_to_nonzero(data, []);
    
    % Data properties
    data_shape = size(data);

    data_properties.spacing = pxs_spacing;
    data_properties.shape_before_cropping = shape_before_cropping;
    data_properties.bbox_used_for_cropping = bbox;
    data_properties.shape_after_cropping_and_before_resampling = data_shape(2:end);
    
    % Adjust target spacing if necessary
    if length(target_spacing) < length(data_shape(2:end))
        target_spacing = [original_spacing(1), target_spacing'];
    end
    
    % Compute new shape
    new_shape = compute_new_shape(data_shape(2:end), original_spacing, target_spacing);
    
    % Normalize data
    for c = 1:size(data, 2)
        %%%%%%%%% CHECK %%%%%%%%%%%%
        normalizer = ZScoreNormalization(use_mask_for_norm(c), IntProp.x0);
        data(c, :, :, :) = normalizer.run(data(c, :, :, :), seg(c, :, :, :));
    end
    
    % Resample data to new shape
    old_shape = data_shape(2:end);
    
    data = resample_data_or_seg_to_shape(data, new_shape, original_spacing, target_spacing);
end


%% 
function segmentation_reverted_cropping = export_prediction_from_softmax_onnx(predicted_array, info_preprocess, properties_dict, erosion)
    
    if nargin < 4
        erosion = 0;
    end

    predicted_array = single(predicted_array);
    % Resample to original shape
    current_spacing = info_preprocess.target_spacing;
    % Resampling
    predicted_probabilities= zeros([size(predicted_array,1),size(predicted_array,2), properties_dict.shape_after_cropping_and_before_resampling], 'single');
    s=size(predicted_array);
    for b=1:size(predicted_probabilities,2)
        tmp_p=reshape(predicted_array(:,b,:,:,:),[s(1) s(3:end)]);
        predicted_probabilities(:,b,:,:,:) = resample_data_or_seg_to_shape(tmp_p, properties_dict.shape_after_cropping_and_before_resampling, current_spacing, properties_dict.spacing);
    end
    % Segmentation does not handle multiple regions
    predicted_probabilities = np_sigmoid(predicted_probabilities);

    [~, segmentation] = max(predicted_probabilities, [], 2); %1 2

    % Adjust 'segmentation' shape to be compatible with python implementation
    % Drop first dimension
    %segmentation = permute(segmentation, [1,3,4,2]);

    % Put result in bbox (revert cropping)
    segmentation_reverted_cropping = zeros([s(1),properties_dict.shape_before_cropping], 'uint8');
    slicer = bounding_box_to_slice(properties_dict.bbox_used_for_cropping);
    segmentation_reverted_cropping(:,slicer{1}{1}:slicer{1}{2}-1, slicer{2}{1}:slicer{2}{2}-1, slicer{3}{1}: slicer{3}{2}-1) = segmentation;
    
    % Conversions necessary as Matlab maps argmax differently
    segmentation_reverted_cropping(segmentation_reverted_cropping == 1) = 0;
    segmentation_reverted_cropping(segmentation_reverted_cropping == 2) = 1;
    
    clear segmentation 
    transpose_vector = [1, info_preprocess.transpose_backward'+2];
    % Revert Transpose
    segmentation_reverted_cropping = permute(segmentation_reverted_cropping, transpose_vector);
    
    if erosion > 0
        se = ones(1, erosion, erosion);
        segmentation_reverted_cropping = imerode(segmentation_reverted_cropping, se);
    end

end

%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%% Other Functions
%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%%--------------------------------------------------------------------------------------------------------------------------------------------------%%



%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%% nnUnet_utils_numpy 

function nonzero_mask = create_nonzero_mask(data)
    assert((length(size(data)) == 4) | (length(size(data)) == 3))
    
    data_shape = size(data);
    nonzero_mask = zeros(data_shape(2:end), 'logical');
    for c = 1:data_shape(1)
        
        % Segment and remove the first dimension which will always be 1
        data_seg = data(c,:,:,:);
        sz = size(data_seg);
        t2 = sz~=1;
        dims = [2:length(size(data_seg))];
        t2(dims)=1;
        out = reshape(data_seg,sz(t2)); %// out has first dimension removed
    
        this_mask = out ~= 0;
        nonzero_mask = nonzero_mask | this_mask;
    
        % imfill only takes in 2D image (remove first dim)
        shape = size(nonzero_mask);
        squeezed = reshape(nonzero_mask,[shape(2:end) 1]);
    
        nonzero_mask = imfill(squeezed, 'holes');
        
        % Resupply first dim
        nonzero_mask = reshape(nonzero_mask, [1, size(nonzero_mask)]);
    end
end


function [data, seg, bbox] = crop_to_nonzero(data, seg, nonzero_label)
    % Crop the data and segmentation map to nonzero regions
    if nargin < 3
        nonzero_label = -1;
    end
    
    nonzero_mask = create_nonzero_mask(data);
    bbox = get_bbox_from_mask(nonzero_mask);
    slicer = bounding_box_to_slice(bbox);
    
    % -1 is needed because matlab includes the end index 
    data = data(:, slicer{1}{1}:slicer{1}{2}-1, slicer{2}{1}:slicer{2}{2}-1,slicer{3}{1}:slicer{3}{2}-1);

    if ~isempty(seg)
        seg = seg(:, slicer{1}{1}:slicer{1}{2}-1, slicer{2}{1}:slicer{2}{2}-1,slicer{3}{1}:slicer{3}{2}-1);
    end
    
    nonzero_mask = nonzero_mask(slicer{1}{1}:slicer{1}{2}-1, slicer{2}{1}:slicer{2}{2}-1,slicer{3}{1}:slicer{3}{2}-1);
    % The [None] in python adds a dim in the front
    nonzero_mask = reshape(nonzero_mask, [1, size(nonzero_mask)]);
    
    if ~isempty(seg)
        %%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%
        seg(seg == 0 & ~nonzero_mask) = nonzero_label;
    else
        nonzero_mask = int8(nonzero_mask);
        nonzero_mask(nonzero_mask == 0) = nonzero_label;
        nonzero_mask(nonzero_mask > 0) = 0;
        seg = nonzero_mask;
    end
end


function reshaped = resize_segmentation(segmentation, new_shape, order)
    % RESIZE_SEGMENTATION Resize segmentation map
    % :param segmentation: Input segmentation map
    % :param new_shape: Target shape
    % :param order: Order of interpolation (default: 3)
    % :return: reshaped: Resized segmentation map
    
    if nargin < 3
        order = 3;
    end
    
    tpe = class(segmentation); 
    unique_labels = unique(segmentation);
    assert(ndims(segmentation) == length(new_shape), 'new shape must have same dimensionality as segmentation');
    
    if order == 0
        %%%%%%%%%% CHECK %%%%%%%%%%%
        reshaped = imresize3(segmentation, new_shape, 'nearest');
        reshaped = cast(reshaped, tpe);
    else
        reshaped = zeros(new_shape, 'like', segmentation);

        for i = 1:length(unique_labels)
            c = unique_labels(i);
            mask = segmentation == c;
            reshaped_multihot = imresize3(single(mask), new_shape, 'nearest');
            reshaped(reshaped_multihot >= 0.5) = c;
        end
    end
end

function do_separate_z = get_do_separate_z(spacing, anisotropy_threshold)
    % GET_DO_SEPARATE_Z Determine if spacing requires separate resampling along Z
    % :param spacing: Spacing values
    % :param anisotropy_threshold: Anisotropy threshold (default: 3)
    % :return: do_separate_z: Whether to separate Z resampling
    
    if nargin < 2
        anisotropy_threshold = 3;
    end
    
    do_separate_z = (max(spacing) / min(spacing)) > anisotropy_threshold;
end

function axis = get_lowres_axis(new_spacing)
    % GET_LOWRES_AXIS Get the axis with the lowest resolution
    % :param new_spacing: Spacing values
    % :return: axis: Axis with the lowest resolution
    axis = find((max(new_spacing) ./ new_spacing) == 1);
end

function new_shape = compute_new_shape(old_shape, old_spacing, new_spacing)
    % COMPUTE_NEW_SHAPE Compute the new shape based on spacing changes
    % :param old_shape: Original shape
    % :param old_spacing: Original spacing
    % :param new_spacing: New spacing
    % :return: new_shape: New shape
    
    assert(length(old_spacing) == length(old_shape), 'old_spacing and old_shape must have the same length');
    assert(length(old_shape) == length(new_spacing), 'old_shape and new_spacing must have the same length');
    new_shape = round((old_spacing ./ new_spacing) .* old_shape);
end

function data_reshaped = resample_data_or_seg_to_spacing(data, current_spacing, new_spacing, is_seg, order, order_z, force_separate_z, separate_z_anisotropy_threshold)
    % RESAMPLE_DATA_OR_SEG_TO_SPACING Resample data or segmentation to new spacing
    % :param data: Input data
    % :param current_spacing: Current spacing
    % :param new_spacing: New spacing
    % :param is_seg: Whether the input is a segmentation map (default: false)
    % :param order: Order of interpolation (default: 3)
    % :param order_z: Order of interpolation along Z (default: 0)
    % :param force_separate_z: Force separate resampling along Z (default: [])
    % :param separate_z_anisotropy_threshold: Anisotropy threshold for separate Z resampling (default: 3)
    % :return: data_reshaped: Resampled data
    
    if nargin < 4
        is_seg = false;
    end
    if nargin < 5
        order = 3;
    end
    if nargin < 6
        order_z = 0;
    end
    if nargin < 7
        force_separate_z = [];
    end
    if nargin < 8
        separate_z_anisotropy_threshold = 3;
    end
    
    if ~isempty(force_separate_z)
        do_separate_z = force_separate_z;
        if force_separate_z
            axis = get_lowres_axis(current_spacing);
        else
            axis = [];
        end
    else
        if get_do_separate_z(current_spacing, separate_z_anisotropy_threshold)
            do_separate_z = true;
            axis = get_lowres_axis(current_spacing);
        elseif get_do_separate_z(new_spacing, separate_z_anisotropy_threshold)
            do_separate_z = true;
            axis = get_lowres_axis(new_spacing);
        else
            do_separate_z = false;
            axis = [];
        end
    end
    
    if ~isempty(axis)
        if length(axis) == 3
            do_separate_z = false;
        elseif length(axis) == 2
            do_separate_z = false;
        end
    end
    
    if ~isempty(data)
        assert(ndims(data) == 4, 'data must be c x y z');
    end
    
    shape = size(data);
    new_shape = compute_new_shape(shape(2:end), current_spacing, new_spacing);
    data_reshaped = resample_data_or_seg(data, new_shape, is_seg, axis, order, do_separate_z, order_z);
end

function data_reshaped = resample_data_or_seg_to_shape(data, new_shape, current_spacing, new_spacing, is_seg, order, order_z, force_separate_z, separate_z_anisotropy_threshold)
    % RESAMPLE_DATA_OR_SEG_TO_SHAPE Resample data or segmentation to new shape
    % :param data: Input data
    % :param new_shape: New shape
    % :param current_spacing: Current spacing
    % :param new_spacing: New spacing
    % :param is_seg: Whether the input is a segmentation map (default: false)
    % :param order: Order of interpolation (default: 3)
    % :param order_z: Order of interpolation along Z (default: 0)
    % :param force_separate_z: Force separate resampling along Z (default: 0)
    % :param separate_z_anisotropy_threshold: Anisotropy threshold for separate Z resampling (default: 3)
    % :return: data_reshaped: Resampled data
    
    if nargin < 5
        is_seg = false;
    end
    if nargin < 6
        order = 3;
    end
    if nargin < 7
        order_z = 0;
    end
    if nargin < 8
        force_separate_z = 0;
    end
    if nargin < 9
        separate_z_anisotropy_threshold = 3;
    end
    
    if ~isempty(force_separate_z)
        do_separate_z = force_separate_z;
        if force_separate_z
            axis = get_lowres_axis(current_spacing);
        else
            axis = [];
        end
    else
        if get_do_separate_z(current_spacing, separate_z_anisotropy_threshold)
            do_separate_z = true;
            axis = get_lowres_axis(current_spacing);
        elseif get_do_separate_z(new_spacing, separate_z_anisotropy_threshold)
            do_separate_z = true;
            axis = get_lowres_axis(new_spacing);
        else
            do_separate_z = false;
            axis = [];
        end
    end
    
    if ~isempty(axis)
        if length(axis) == 3
            do_separate_z = false;
        elseif length(axis) == 2
            do_separate_z = false;
        end
    end
    
    if ~isempty(data)
        assert(ndims(data) == 4, 'data must be c x y z');
    end
    
    data_reshaped = resample_data_or_seg(data, new_shape, is_seg, axis, order, do_separate_z, order_z);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reshaped_final_data = resample_data_or_seg(data, new_shape, is_seg, axis, order, do_separate_z, order_z)
    % RESAMPLE_DATA_OR_SEG Resample data or segmentation
    % :param data: Input data
    % :param new_shape: New shape
    % :param is_seg: Whether the input is a segmentation map (default: false)
    % :param axis: Axis for separate Z resampling (default: [])
    % :param order: Order of interpolation (default: 3)
    % :param do_separate_z: Whether to separate Z resampling (default: false)
    % :param order_z: Order of interpolation along Z (default: 0)
    % :return: reshaped_final_data: Resampled data
    
    if nargin < 3
        is_seg = false;
    end
    if nargin < 4
        axis = [];
    end
    if nargin < 5
        order = 3;
    end
    if nargin < 6
        do_separate_z = false;
    end
    if nargin < 7
        order_z = 0;
    end
    
    assert(ndims(data) == 4, 'data must be (c, x, y, z)');
    assert(length(new_shape) == ndims(data) - 1, 'new shape must match data dimensions minus one');
    
    if is_seg
        resize_fn = @resize_segmentation;
    else
        resize_fn = @imresize3;
        kwargs.mode = 'edge';
        kwargs.anti_aliasing = false;
    end
    
    dtype_data = class(data);
    data_slice = data(1,:,:,:); 
    data_slice = permute(data_slice, [linspace(2, ndims(data), ndims(data)-1), 1]);
    shape = [size(data_slice)];
    new_shape = [new_shape];
    
    if any(shape ~= new_shape)
        data = double(data);
        if do_separate_z
            assert(length(axis) == 1, 'only one anisotropic axis supported');
            axis = axis(1);
            if axis == 1
                new_shape_2d = new_shape(2:end);
            elseif axis == 2
                new_shape_2d = new_shape([1, 3]);
            else
                new_shape_2d = new_shape(1:end-1);
            end
            
            reshaped_final_data = [];
            for c = 1:size(data, 1)
                reshaped_data = [];
                for slice_id = 1:shape(axis)
                    if axis == 1
                        if is_seg 
                            reshaped_slice = resize_fn(squeeze(data(c, slice_id, :, :)), new_shape_2d, order);
                        else 
                            reshaped_slice = resize_fn(squeeze(data(c, slice_id, :, :)), new_shape_2d, 'nearest');
                        end
                    
                    elseif axis == 2
                        if is_seg
                            reshaped_slice = resize_fn(squeeze(data(c, :, slice_id, :)), new_shape_2d, order);
                        else
                            reshaped_slice = resize_fn(squeeze(data(c, :, slice_id, :)), new_shape_2d, 'nearest');
                        end
                    else
                        if is_seg
                            reshaped_slice = resize_fn(squeeze(data(c, :, :, slice_id)), new_shape_2d, order);
                        else
                            reshaped_slice = resize_fn(squeeze(data(c, :, :, slice_id)), new_shape_2d, order,  'nearest');
                        end
                    end
                    reshaped_data = cat(axis, reshaped_data, reshaped_slice);
                end
                
                if shape(axis) ~= new_shape(axis)
                    % Conversion to assign multiple variables at once
                    new_shape_cell = num2cell(new_shape);
                    [rows, cols, dim] = deal(new_shape_cell{:});

                    reshaped_data_shape = size(reshaped_data);
                    reshaped_data_shape_cell = num2cell(reshaped_data_shape);

                    [orig_rows, orig_cols, orig_dim] = deal(reshaped_data_shape_cell{:});
                    
                    row_scale = orig_rows / rows;
                    col_scale = orig_cols / cols;
                    dim_scale = orig_dim / dim;
                    
                    [map_rows, map_cols, map_dims] = ndgrid(1:rows, 1:cols, 1:dims);
                    map_rows = row_scale * (map_rows + 0.5) - 0.5;
                    map_cols = col_scale * (map_cols + 0.5) - 0.5;
                    map_dims = dim_scale * (map_dims + 0.5) - 0.5;

                    coord_map = cat(4, map_rows, map_cols, map_dims);
                    coord_map = permute(coord_map, [4,1,2,3]); 

                    if ~is_seg || order_z == 0
                        %%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%
                        reshaped_final_data = cat(1, reshaped_final_data, reshape(imresize3(reshaped_data, coord_map, 'nearest'), [1, size(reshaped_data)]));
                    else
                        unique_labels = unique(reshaped_data(:));
                        reshaped = zeros(new_shape, 'like', dtype_data);
                        
                        for i = 1:length(unique_labels)
                            cl = unique_labels(i);
                            reshaped_multihot = round(imresize3(single(reshaped_data == cl), coord_map, 'nearest'));
                            reshaped(reshaped_multihot > 0.5) = cl;
                        end
                        reshaped_final_data = cat(1, reshaped_final_data, reshape(reshaped, [1, size(reshaped)]));
                    end
                else
                    reshaped_final_data = cat(1, reshaped_final_data, reshape(reshaped_data, [1, size(reshaped_data)]));
                end
            end
        else
            reshaped = [];
            for c = 1:size(data, 1)
                reshaped = cat(1, reshaped, reshape(resize_fn(reshape(squeeze(data(c, :, :, :)),shape), new_shape, 'linear'), [1, new_shape]));
            end
            reshaped_final_data = reshaped;
        end
        
        reshaped_final_data = cast(reshaped_final_data, dtype_data);
    else
        reshaped_final_data = data;
    end
end
  

%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%% acvl_utils_bounding_boxes 

function bbox_utils()
    bbox = [32, 64; 21, 46];
    bbox_padded = pad_bbox(bbox, 3);
    slicer = bounding_box_to_slice(bbox_padded);
end

function bounding_box = pad_bbox(bounding_box, pad_amount, array_shape)
    if nargin < 3
        array_shape = [];
    end
    
    if isstruct(bounding_box)
        bounding_box = struct2cell(bounding_box);
        bounding_box = cell2mat(bounding_box);
    elseif isnumeric(bounding_box)
        bounding_box = double(bounding_box);
    else
        error('bounding_box must be a numeric array or struct');
    end

    if isscalar(pad_amount)
        pad_amount = repmat(pad_amount, size(bounding_box, 1), 1);
    else
        pad_amount = double(pad_amount);
    end

    for i = 1:size(bounding_box, 1)
        new_values = [max(0, bounding_box(i, 1) - pad_amount(i)), bounding_box(i, 2) + pad_amount(i)];
        if ~isempty(array_shape)
            new_values(2) = min(array_shape(i), new_values(2));
        end
        bounding_box(i, :) = new_values;
    end
end

function proper_bbox = regionprops_bbox_to_proper_bbox(regionprops_bbox)
    dim = numel(regionprops_bbox) / 2;
    proper_bbox = zeros(dim, 2);
    for i = 1:dim
        proper_bbox(i, :) = [regionprops_bbox(i), regionprops_bbox(i + dim)];
    end
end

function slicer = bounding_box_to_slice(bounding_box)
    slicer = cell(1, size(bounding_box, 2));
    for i = 1:size(bounding_box, 2)
        slicer{i} = bounding_box{i};
    end
end

function cropped_array = crop_to_bbox(array, bounding_box)
    assert(size(bounding_box, 1) == ndims(array), ...
        sprintf('Dimensionality of bbox and array do not match. bbox has length %d while array has dimension %d', ...
        size(bounding_box, 1), ndims(array)));
    slicer = bounding_box_to_slice(bounding_box);
    cropped_array = array(slicer{:});
end

function bbox = get_bbox_from_mask(nonzero_mask)
    [Z, X, Y] = size(nonzero_mask);
    values = {0, Z, 0, X, 0, Y};
    [minzidx, maxzidx, minxidx, maxxidx, minyidx, maxyidx] = values{:};
    
    zidx = linspace(1,Z,Z);
    
    for z = zidx
        mask = nonzero_mask(z,:,:);
    
        % Drop the first dim
        shape = size(mask);
        mask = reshape(mask, [shape(2:end), 1]);
        
        if any(mask(:)) 
            minzidx = z;
            break 
        end
    end
    
    for z = linspace(Z,1,Z)
        mask = nonzero_mask(z,:,:);
    
        % Drop the first dim
        shape = size(mask);
        mask = reshape(mask, [shape(2:end), 1]);
    
        if any(mask(:)) 
            maxzidx = z + 1;
        end
    end
    
    xidx = linspace(1,X,X);
    for x = xidx
        mask = nonzero_mask(:,x,:);
    
        % Drop the first dim
        shape = size(mask);
        mask = reshape(mask, [shape(2:end), 1]);
       
        
        if any(mask(:)) 
            minxidx = x;
            break 
        end
    end
    
    for x = linspace(X,1,X)
        mask = nonzero_mask(:,x,:);
    
        % Drop the first dim
        shape = size(mask);
        mask = reshape(mask, [shape(2:end), 1]);
       
        if any(mask(:)) 
            maxxidx = x + 1;
            break
        end
    end
    
    yidx = linspace(1,Y,Y);
    for y = yidx
        mask = nonzero_mask(:,:,y);
        
        if any(mask(:)) 
            minyidx = y;
            break 
        end
    end
    
    for y = linspace(Y,1,Y)
        mask = nonzero_mask(:,:,y);  
    
        if any(mask(:)) 
            maxyidx = y + 1;
            break
        end
    end
    bbox = {{minzidx, maxzidx}, {minxidx, maxxidx}, {minyidx, maxyidx}};
                  
end

function bbox = get_bbox_from_mask_npwhere(mask)
    [Z, X, Y] = size(mask);
    [z, x, y] = ind2sub(size(mask), find(mask));
    mins = [min(z), min(x), min(y)];
    maxs = [max(z), max(x), max(y)] + 1;
    bbox = [mins', maxs'];
end


%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%% Sliding Window 

% edited
function gaussian_importance_map = compute_gaussian(tile_size, sigma_scale)
    if nargin < 2
        sigma_scale = 1/8;
    end
    tmp = zeros(tile_size');
    center_coords = floor(tile_size/ 2) + 1;
    sigmas = tile_size * sigma_scale;
    tmp(center_coords(1), center_coords(2)) = 1;

    gaussian_importance_map = fspecial('gaussian',[tile_size(1)+1, tile_size(2)+1] ,sigmas(1));
    gaussian_importance_map = gaussian_importance_map(1:end-1, 1:end-1);
    gaussian_importance_map = gaussian_importance_map / max(gaussian_importance_map(:));
    gaussian_importance_map = cast(gaussian_importance_map, 'half');
    gaussian_importance_map(gaussian_importance_map == 0) = min(gaussian_importance_map(gaussian_importance_map ~= 0));
end


function steps = compute_steps_for_sliding_window(image_size, tile_size, tile_step_size)
    %image_size = image_size(2:end);
    % disp(image_size)
    assert(all(image_size >= tile_size), 'image size must be as large or larger than tile_size');
    assert(tile_step_size > 0 & tile_step_size <= 1, 'step_size must be larger than 0 and smaller or equal to 1');

    target_step_sizes_in_voxels = tile_size * tile_step_size;
    num_steps = ceil((image_size - tile_size) ./ target_step_sizes_in_voxels) + 1;

    steps = cell(1, length(tile_size));
    for dim = 1:length(tile_size)
        max_step_value = image_size(dim) - tile_size(dim);
        if num_steps(dim) > 1
            actual_step_size = max_step_value / (num_steps(dim) - 1);
        else
            actual_step_size = 99999999999;
        end
        steps{dim} = round(actual_step_size * (0:(num_steps(dim) - 1)));
    end
end

function slicers = get_sliding_window_generator(image_size, tile_size, tile_step_size, verbose)
    if nargin < 4
        verbose = false;
    end

    % Transpose to make operations consistent
    tile_size = tile_size'; 
    tile_step_size = tile_step_size';

    if length(tile_size) < length(image_size)
        assert(length(tile_size) == length(image_size) - 1, 'if tile_size has less entries than image_size, len(tile_size) must be one shorter than len(image_size)');
                
        steps = compute_steps_for_sliding_window(image_size(2:end), tile_size, tile_step_size);
    
        if verbose
            disp(['n_steps ', num2str(image_size(1) * prod(cellfun(@length, steps))), ', image size is ', mat2str(image_size), ', tile_size ', mat2str(tile_size), ', tile_step_size ', num2str(tile_step_size), ' steps: ', mat2str(cell2mat(steps(:)'))]);
        end
        slicers = cell(1, image_size(1) * prod(cellfun(@length, steps)));
        idx = 1;
        for d = 1:image_size(1)
            for sx = steps{1}
                for sy = steps{2}
                    slicers{idx} = {':', d, sx + (1:sx+tile_size(1)), sy + (1:sy+tile_size(2))};
                    idx = idx + 1;
                end
            end
        end
    
    else 
        steps = compute_steps_for_sliding_window(image_size, tile_size, tile_step_size);

        if verbose
            disp(['n_steps ', num2str(prod(cellfun(@length, steps))), ', image size is ', mat2str(image_size), ', tile_size ', mat2str(tile_size), ', tile_step_size ', num2str(tile_step_size), ' steps: ', mat2str(cell2mat(steps(:)'))]);
        end
        slicers = cell(1, prod(cellfun(@length, steps)));
        idx = 1;
        for sx = steps{1}
            for sy = steps{2}
                for sz = steps{3}
                    slicers{idx} = {':', sx + (1:sx+tile_size(1)), sy + (1:sy+tile_size(2)), sz + (1:sz+tile_size(3))};
                    idx = idx + 1;
                end
            end
        end
    end
end

function prediction = maybe_mirror_and_predict(network, x, mirror_axes)
    if nargin < 3
        mirror_axes = [];
    end
    %% X : BCSS
    %PD
    %input_name = network.InputNames{1};
    %x = cast(x,'single');
    %ort_inputs = containers.Map(input_name, x);
    %prediction = network.predict(squeeze(ort_inputs(input_name)));
    % To achieve consistent dimensions with python implementation
    %prediction = reshape(prediction, [1, size(prediction)]);
    %prediction = permute(prediction, [1,4,2,3]);
    x_mat=permute(x,[4,3,2,1]);
    input_x=dlarray(single(x_mat),"SSCB");
    if (canUseGPU) 
    input_x = gpuArray(input_x);
    end
    dlprediction = network.predict(input_x);
    prediction =extractdata(dlprediction); %%SSCB
    prediction = permute(prediction, [4,3,2,1]);
    

    if ~isempty(mirror_axes)
        assert(max(mirror_axes) <= ndims(x) - 3, 'mirror_axes does not match the dimension of the input!');
    
        num_predictions = 2 ^ length(mirror_axes);
    
        if ismember(0, mirror_axes)
            %PD
            %ort_inputs = containers.Map(input_name, flip(x, 3));
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,4,2,3]);
            x_mat=permute(flip(x,3),[4,3,2,1]);
            input_x=dlarray(single(x_mat),"SSCB");
            if (canUseGPU) 
            input_x = gpuArray(input_x);
            end
            dlprediction = network.predict(input_x);
            predict =extractdata(dlprediction); %%SSCB
            predict = permute(predict, [4,3,2,1]);
    
            prediction = prediction + flip(predict, 3);
        end
    
        if ismember(1, mirror_axes)
            %ort_inputs = containers.Map(input_name, flip(x, 4));
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,4,2,3]);
            x_mat=permute(flip(x,4),[4,3,2,1]);
            input_x=dlarray(single(x_mat),"SSCB");
            if (canUseGPU) 
            input_x = gpuArray(input_x);
            end
            dlprediction = network.predict(input_x);
            predict =extractdata(dlprediction); %%SSCB
            predict = permute(predict, [4,3,2,1]);
            prediction = prediction + flip(predict, 4);
            
        end
        % NEED TO CHECK % 
        if ismember(2, mirror_axes)
            %ort_inputs = containers.Map(input_name, flip(x, 5));
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,5,2,3,4]);
            %x_mat=permute(flip(x,5),[4,3,2,1]);
            %input_x=dlarray(single(x_mat),"SSCB");
            %dlprediction = network.predict(input_x);
            %predict =extractdata(dlprediction); %%SSCB
            %predict = permute(predict, [4,3,2,1]);
            %prediction = prediction + flip(predict, 5);
            disp('Not implemented and tested')
        end
    
        if ismember(0, mirror_axes) && ismember(1, mirror_axes)
            % Matlab doesnt allow inputting flip array
            %input_x = flip(x,3);
            %input_x = flip(input_x,4);
            %ort_inputs = containers.Map(input_name, input_x);
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,4,2,3]);
            
            %predict = flip(predict, 3);
            %predict = flip(predict, 4);
            x_mat=permute(flip(flip(x,3),4),[4,3,2,1]);
            input_x=dlarray(single(x_mat),"SSCB");
            if (canUseGPU) 
            input_x = gpuArray(input_x);
            end
            dlprediction = network.predict(input_x);
            predict =extractdata(dlprediction); %%SSCB
            predict = permute(predict, [4,3,2,1]);

            prediction = prediction + flip(flip(predict,3),4);
        end
    
        if ismember(0, mirror_axes) && ismember(2, mirror_axes)
            disp('Not implemented and tested')

            %input_x = flip(x,3);
            %input_x = flip(input_x, 5);
    
            %ort_inputs = containers.Map(input_name, input_x);
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,5,2,3,4]);
    
            %predict = flip(predict, 3);
            %predict = flip(predict, 5);
            %prediction = prediction + predict;
        end
    
        if ismember(1, mirror_axes) && ismember(2, mirror_axes)
            disp('Not implemented and tested')
            %input_x = flip(x,4);
            %input_x = flip(input_x, 5);
    
            %ort_inputs = containers.Map(input_name, input_x);
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,5,2,3,4]);
    
            %predict = flip(predict, 4);
            %predict = flip(predict, 5);
            %prediction = prediction + predict;
        end
    
        if ismember(0, mirror_axes) && ismember(1, mirror_axes) && ismember(2, mirror_axes)
            disp('Not implemented and tested')
            %input_x = flip(x,3);
            %input_x = flip(input_x, 4);
            %input_x = flip(input_x, 5);
    
            %ort_inputs = containers.Map(input_name, input_x);
            %predict = network.predict(squeeze(ort_inputs(input_name)));
            %predict = reshape(predict, [1, size(predict)]);
            %predict = permute(predict, [1,5,2,3,4]);
    
            %predict = flip(predict, 3);
            %predict = flip(predict, 4);
            %predict = flip(predict, 5);
            %prediction = prediction + predict;
        end

        prediction = prediction / num_predictions;
        if (canUseGPU) 
        prediction=gather(prediction);
        end
    end
end


function predicted_logits = predict_sliding_window_return_logits(networks, input_image, num_segmentation_heads, tile_size, varargin)
    % Default parameter values
    mirror_axes = [];
    tile_step_size = 0.5;
    use_gaussian = true;
    precomputed_gaussian = [];
    verbose = false;

    % Parse optional parameters
    if ~isempty(varargin)
        % Ensure that varargin contains pairs of parameter names and values
        if mod(length(varargin), 2) ~= 0
            error('Optional parameters should be in name-value pairs.');
        end

        for i = 1:2:length(varargin)
            param_name = varargin{i};
            param_value = varargin{i+1};

            if ~ischar(param_name)
                error('Parameter name must be a character vector.');
            end

            switch param_name
                case 'mirror_axes'
                    mirror_axes = param_value;
                case 'tile_step_size'
                    tile_step_size = param_value;
                case 'use_gaussian'
                    use_gaussian = param_value;
                case 'precomputed_gaussian'
                    precomputed_gaussian = param_value;
                case 'verbose'
                    verbose = param_value;
                otherwise
                    error(['Unknown parameter name: ', param_name]);
            end
        end
    end

    assert(ndims(input_image) == 4, 'input_image must be a 4D array');

    % Pad input_image if it is smaller than tile_size
    [data, slicer_revert_padding] = pad_nd_image(input_image, tile_size, 'constant', struct('constant_values', 0), true, []);
    
    if use_gaussian
        if isempty(precomputed_gaussian)
            gaussian = compute_gaussian(tile_size, 1/8);
        else
            gaussian = precomputed_gaussian;
        end
        mn = min(gaussian(:));
        if mn == 0
            gaussian = max(gaussian, mn);
        end
    else
        gaussian = ones(tile_size);
    end

    data_shape = size(data);
    slicers = get_sliding_window_generator(data_shape(2:end), tile_size, tile_step_size, verbose); 
      
    predicted_logits = zeros([data_shape(1),num_segmentation_heads, data_shape(2:end)], 'half');
    n_predictions = zeros(data_shape, 'half');

    % Adjust Gaussian, predicted_logits, and n_predictions dimensions for operations below
    gaussian = repmat(reshape(gaussian, [1,1,size(gaussian)]),[data_shape(1),1,1,1]); %Does not work for 5D data Be careful

    for sl = slicers
        workon = data(slicers{1}{1}, slicers{1}{2}, slicers{1}{3}, slicers{1}{4});
        % Second dim is always dropped
        %workon = workon(:,1,:,:);
    
        % Single dimension added to the front (use the singleton dropped dim)
        %workon = permute(workon, [2,1,3,4]); PD

        for i = 1:length(networks)
            tic
            prediction = maybe_mirror_and_predict(networks{i}, workon, mirror_axes);
            toc
            if use_gaussian
                prediction=prediction.*gaussian;
            end
            predicted_logits(slicers{1}{1},slicers{1}{1},slicers{1}{2}, slicers{1}{3}, slicers{1}{4}) = squeeze(predicted_logits(slicers{1}{1},slicers{1}{1},slicers{1}{2}, slicers{1}{3}, slicers{1}{4})) + squeeze(prediction);
            n_predictions(slicers{1}{1}, slicers{1}{2}, slicers{1}{3}, slicers{1}{4})=squeeze(n_predictions(slicers{1}{1}, slicers{1}{2}, slicers{1}{3}, slicers{1}{4}))+squeeze(gaussian);
        end
    end

    for i = 1:size(predicted_logits,2)
        predicted_logits(:,i,:,:,:) = squeeze(predicted_logits(:,i,:,:,:)) ./ squeeze(n_predictions);
    end
    
    slicer_revert_padding_bis = slicer_revert_padding;
    slicer_revert_padding_bis{1} = ':';
    slicer_revert_padding_bis{2} = ':'; 
    for i =1:size(slicer_revert_padding,2)-1
        slicer_revert_padding_bis{2+i} =slicer_revert_padding{i+1};
    end
    predicted_logits = predicted_logits(slicer_revert_padding_bis{1}, slicer_revert_padding_bis{2}, slicer_revert_padding_bis{3}, slicer_revert_padding_bis{4},slicer_revert_padding_bis{5});
end


%%--------------------------------------------------------------------------------------------------------------------------------------------------%%
%% acvl utils numpy, edited

function [res, slicer] = pad_nd_image(image, new_shape, mode, kwargs, return_slicer, shape_must_be_divisible_by)
     
     % Handle default values for the parameters
    if nargin < 2, new_shape = []; end
    if nargin < 3, mode = 'constant'; end
    if nargin < 4, kwargs = struct(); end
    if nargin < 5, return_slicer = false; end
    if nargin < 6, shape_must_be_divisible_by = []; end
    
    old_shape = size(image);

    % Handle shape_must_be_divisible_by
    if ~isempty(shape_must_be_divisible_by)
        assert((isnumeric(shape_must_be_divisible_by) && all(mod(shape_must_be_divisible_by(:), 1) == 0)) || ...
        (iscell(shape_must_be_divisible_by) && all(cellfun(@(x) isnumeric(x) && all(mod(x(:), 1) == 0), shape_must_be_divisible_by))), ...
        'shape_must_be_divisible_by must be an integer, a list of integers, a tuple of integers, or a numeric array of integers.');
        if length(shape_must_be_divisible_by) == 1 && floor(shape_must_be_divisible_by) == shape_must_be_divisible_by
            shape_must_be_divisible_by = repmat(shape_must_be_divisible_by, 1, ndims(image));
        else 
            if length(shape_must_be_divisible_by) < ndims(image)
                shape_must_be_divisible_by = [ones(1, ndims(image) - length(shape_must_be_divisible_by)), shape_must_be_divisible_by];
            end
        end
    end
    
    % Handle new_shape
    if isempty(new_shape)
        assert(~isempty(shape_must_be_divisible_by), 'shape_must_be_divisible_by must be provided if new_shape is not.');
        new_shape = size(image);
    end

    if length(new_shape) < ndims(image)
        image_shape = size(image);
        new_shape = [image_shape(1: ndims(image) - length(new_shape)), new_shape'];
    end

    new_shape = arrayfun(@(i) max(new_shape(i), old_shape(i)), 1:length(new_shape));
    
    % Adjust new_shape based on shape_must_be_divisible_by
    if ~isempty(shape_must_be_divisible_by)
        if length(shape_must_be_divisible_by) == 1
            shape_must_be_divisible_by = repmat(shape_must_be_divisible_by, 1, length(new_shape));
        end
        if length(shape_must_be_divisible_by) < length(new_shape)
            shape_must_be_divisible_by = [ones(1, length(new_shape) - length(shape_must_be_divisible_by)), shape_must_be_divisible_by];
        end
    

        for i = 1:length(new_shape)
            if mod(new_shape(i), shape_must_be_divisible_by(i)) == 0
                new_shape(i) = new_shape(i) - shape_must_be_divisible_by(i);
            end
            new_shape(i) = new_shape(i) + shape_must_be_divisible_by(i) - mod(new_shape(i), shape_must_be_divisible_by(i));
        end
    end
    
    difference = new_shape - old_shape;
    pad_below = floor(difference / 2);
    pad_above = ceil(difference / 2);
    pad_list = [pad_below(:), pad_above(:)];
    
    
    if ~all(pad_below == 0 & pad_above == 0)
        res = padarray(image, pad_below, 0, 'pre');
        res = padarray(res, pad_above, 0, 'post');
    else
        res = image;
    end
    
    if ~return_slicer
        slicer = [];
    else
        pad_list(:, 2) = size(res)' - pad_list(:, 2);
        slicer = arrayfun(@(i) pad_list(i, 1)+1:pad_list(i, 2), 1:size(pad_list, 1), 'UniformOutput', false);
        %slicer = cell2mat(slicer);
    end
end
