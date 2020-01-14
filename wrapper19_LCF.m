%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program for fitting ceria XANES data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1;
ratio_19 = [];
error_19 = [];
sample_19 = [];
fit = [];

%% Import the data
Ce3plus_ref = importdata('Ce_NO3_3_merge.txt');
Ce4plus_ref = importdata('CeO2_3_merge.txt');
for i = 1:length(file)
[samplex sampley] = textread(['G:\RT_vac\',file(i).name],'%f %f');
%% Plot the data
if(verbose)
    figure;
    plot(Ce3plus_ref.data(:,1), Ce3plus_ref.data(:, 2), 'r-', ...
         Ce4plus_ref.data(:,1), Ce4plus_ref.data(:, 2), 'b-', ...
         samplex, sampley, 'k--');
    xlabel('Energy (eV)');
    ylabel('\mu / absorption coefficient');
    legend('3+', '4+', 'sample');
    title('Raw reference and sample data');
end
%% Calculate the start and stop values for the interpolation
start_value = ceil(max(min([Ce3plus_ref.data(:,1); Ce3plus_ref.data(:,1); samplex])));
end_value = floor(min(max([Ce3plus_ref.data(:,1); Ce3plus_ref.data(:,1); sampley])));
start_value = 5600;
end_value = 5900;
%% Interpolate the data to get everything on the same X-axis
new_range = start_value:0.25:end_value;
Ce3p = interp1(Ce3plus_ref.data(:, 1), Ce3plus_ref.data(:, 2), new_range);
if(verbose) 
    figure;
    plot(new_range, Ce3p, 'o', Ce3plus_ref.data(:, 1), Ce3plus_ref.data(:, 2), 'k')
    title('Interpolated data')
end

Ce4p = interp1(Ce4plus_ref.data(:, 1), Ce4plus_ref.data(:, 2), new_range);
smpld = interp1(samplex, sampley, new_range);

%% Use the linear combination function to fit the sample data to reference

[p_optimal, sumsq, resid, eflag, op, ~, jacob] = lsqcurvefit(@(x,xdata) lincomb(x, Ce3p, Ce4p), [0.1, 0.9, 0], ...
                new_range, smpld);
x1 = p_optimal(1);
x2 = p_optimal(2);
fit_curve = lincomb(p_optimal, Ce3p, Ce4p);

% Calculate the confidence intervals
ci_95 = nlparci(p_optimal, resid, 'jacobian', jacob);
error = (ci_95(:, 2) - p_optimal');
error_pct = error ./ abs(p_optimal)';

ratio_19 = [ratio_19; x1 x2];
error_19 = [error_19; error(1) error(2)];
fit = [fit fit_curve'];
%%
figure;
plot(new_range, fit_curve, 'b--', ...
     samplex, sampley, 'ok');
xlabel('Energy (eV)');
ylabel('\mu / absorption coefficient');
legend('fit', 'sample');
title('Fitted data');

%% Plot confidence in error
figure;
hold on;
errorbar(0.45, x1, error(1), error(1), 'mx');
errorbar(0.45, x2, error(2), error(2), 'bx');
legend('Ce^{3+}', 'Ce^{4+}');
title('Fit parameters');
xlabel('composition');
ylabel('parameter value');
i = i + 1;
end
fit_19 = [new_range' fit];
save(filename,ratio,curve )