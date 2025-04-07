classdef ImageNormalization < handle
    properties
        use_mask_for_norm
        intensityproperties
        target_dtype
    end
    
    properties (Abstract, Constant)
        leaves_pixels_outside_mask_at_zero_if_use_mask_for_norm_is_true
    end
    
    methods
        function obj = ImageNormalization(use_mask_for_norm, intensityproperties, target_dtype)
            if nargin < 1
                use_mask_for_norm = [];
            end
            if nargin < 2
                intensityproperties = struct();
            end
            if nargin < 3
                target_dtype = 'single';
            end
            
            assert(islogical(use_mask_for_norm) || isempty(use_mask_for_norm), ...
                'use_mask_for_norm must be a boolean or empty');
            obj.use_mask_for_norm = use_mask_for_norm;

            assert(isstruct(intensityproperties), 'intensityproperties must be a struct');
            obj.intensityproperties = intensityproperties;
            obj.target_dtype = target_dtype;
        end
        
        function result = run(obj, image, seg)
            % Abstract method
            error('run method must be implemented by subclasses');
        end
    end
end