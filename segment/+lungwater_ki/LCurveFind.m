function lambda = LCurveFind(image, lambda_range)


if nargin < 2
   lambda_range=linspace(0.5,2000); %lambda search range
else
    lambda_range=linspace(lambda_range(1),lambda_range(2)); %lambda search range
end
image=double(image);

ind = image>0; %values less than zero will be excluded from the fitting process
b = image(image>0); 


% lambda_range=linspace(0.5,2000); %lambda search range



for i = 1:length(lambda_range)
    %Fitting the surface with a specific smoothing parameter
    [x,A,T] = lungwater_ki.tikReg2D(image,lambda_range(i));

    %calculating the errors with that parameter
    res_norm(i) = norm(A*x(:)-b,'fro');
    solution_norm(i) = norm(T*x(:),'fro');
end

res_norm_log = log(res_norm);
solution_norm_log = log(solution_norm);

x_grid = 0.5:0.25:50000;
%interpolate norms
res_norm_log= spline(lambda_range,res_norm_log,x_grid);
solution_norm_log = spline(lambda_range,solution_norm_log,x_grid);

%calculating maximum curvature, derivatives
xL1 = gradient(res_norm_log);
yL1 = gradient(solution_norm_log);

xL2 = del2(res_norm_log);
yL2 = del2(solution_norm_log);

k = (xL2.*yL1-xL1.*yL2)./(xL1.^2+yL1.^2).^1.5; %curvature equations
[~,ind] = min(k);
lambda = x_grid(ind); %optimized lambda at max curvature

end







