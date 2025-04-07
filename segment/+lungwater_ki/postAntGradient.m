function [slope,sliceProfilePosteriorAnterior]= postAntGradient(no, nbr_of_regions, doPlot, filename)

global SET

if nargin<2
    nbr_of_regions=10;
end

if nargin <3
    doPlot=1;
end

if nargin<4
    filename='';
end
tf=1;

LungMask=squeeze(SET(no).LungWater.LungMask(:,:,tf,:));
slices=sort(unique(find(squeeze(sum(sum(LungMask)))>0)))';  %slices with mask

LungMask(LungMask==0)=NaN;
LWD_map=squeeze(SET(no).LungWater.LWD_map(:,:,tf,:));

%divide into regions
slices_per_group=floor(length(slices)/nbr_of_regions);
if slices_per_group==0
    slices_per_group=1;
end

startInd=1;
for i=1:nbr_of_regions
    if i==nbr_of_regions
        slice_groups(i, 1:2) = [slices(startInd) slices(end)];
    else
        slice_groups(i,1:2) = [slices(startInd) slices(startInd+slices_per_group-1)];
        startInd=i*floor(length(slices)/nbr_of_regions)+1;
    end
end

%mean LWD per region
for i=1:nbr_of_regions
    sliceProfilePosteriorAnterior(i,1)=nanmean(nanmean(nanmean(LWD_map(:,:,slice_groups(i,1):slice_groups(i,2)).*LungMask(:,:,slice_groups(i,1):slice_groups(i,2)))));
end

sliceProfilePosteriorAnterior=flip(sliceProfilePosteriorAnterior); %flip to get positive slope value

%first and second degree polynomial fit
x=linspace(slices_per_group,length(slices)*SET(no).ResolutionY, nbr_of_regions);
p1=polyfit(x, sliceProfilePosteriorAnterior,1);
p2=polyfit(x, sliceProfilePosteriorAnterior,2);
slope=p1(1); %slope from the linear function that is added as a lung water parameter

SET(no).LungWater.AnteriorPosteriorLinearGradient=slope;

%getting R and p value from linear function 
fittedValues = polyval(p1, x);
[rho_linear, pval_linear] = corr(fittedValues(:), sliceProfilePosteriorAnterior(:));
  
%exponential fit with offset
sliceProfilePosteriorAnteriorNormalized = sliceProfilePosteriorAnterior./sliceProfilePosteriorAnterior(10);
customEqn = @(a, b, c, x) b * exp(a * x) + c;
[f1_offset, gof, output_s] = fit(x(:), sliceProfilePosteriorAnteriorNormalized(:), customEqn, 'StartPoint', [0, 0, 0]);


% Get the coefficient values (a, b, c)
coeff_values_f1 = coeffvalues(f1_offset);


% Extract parameter a
a_f1 = coeff_values_f1(1);
b_f1 = coeff_values_f1(2);
c_f1 = coeff_values_f1(3);


SET(no).LungWater.sliceProfilePosteriorAnterior=sliceProfilePosteriorAnterior;
SET(no).LungWater.sliceProfilePosteriorAnteriorNormalized=sliceProfilePosteriorAnteriorNormalized;

SET(no).LungWater.ExponentialRMS1 = gof.rmse;
SET(no).LungWater.ExponentialRsquare = gof.rsquare;
SET(no).LungWater.Exponentiala1 = a_f1;
SET(no).LungWater.Exponentialb1 = b_f1;
SET(no).LungWater.Exponentialc1 = c_f1;


if doPlot
    figure; clf; hold on;
    %plot(x, sliceProfilePosteriorAnteriorNormalized, 'o')
    plot(x, sliceProfilePosteriorAnteriorNormalized, 'o')
    %plot(x, sliceProfilePosteriorAnteriorNormalizedPercentage, 'o')
    hold on;
    %plot(x, p1(1)*x+p1(2), 'r-')
    %plot(x, p2(1)*x.^2+p2(2)*x+p2(3), 'm-')
    xlabel('Lung width (mm)')
    ylabel('Mean LWD (%)')

    %plot(f1, 'b-');
    %plot(f2, 'g-');
    plot(f1_offset, 'b-')
    %plot(f2_offset, 'r-')

    equation1 = sprintf('f(x) = %.2f * exp(%.2f * x) + %.2f', b_f1, a_f1, c_f1);
    %ecquation2 = sprintf('f(x) = %.2f * exp(%.2f * x) + %.2f', b_f2, a_f2, c_f2);

    legend('slice/posteriorSlice', equation1, 'Location', 'Best')
    
    title(filename);

    %legend(equation, 'Location', 'Best')

    %legend('anterior-posterior bins','f(x) = a * x + b','f(x) = a * x^2 + b * x + c','f(x) = a * exp(b * x)', 'f_2(x) = a * exp(b * x) + c * exp(d * x)', 'Location', 'Best')
end
