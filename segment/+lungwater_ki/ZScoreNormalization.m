classdef ZScoreNormalization < ImageNormalization
    properties (Constant)
        leaves_pixels_outside_mask_at_zero_if_use_mask_for_norm_is_true = true;
    end
    
    methods
        function result = run(obj, image, seg)
            if nargin < 3
                seg = [];
            end
            
            image = cast(image, obj.target_dtype);
            
            if ~isempty(obj.use_mask_for_norm) && obj.use_mask_for_norm
                mask = seg >= 0;
                mean_val = mean(image(mask));
                std_val = std(image(mask));
                image(mask) = (image(mask) - mean_val) / max(std_val, 1e-8);
            else
                mean_val = mean(image(:));
                std_val = std(image(:));
                image = (image - mean_val) / max(std_val, 1e-8);
            end
            
            result = image;
        end
    end
end