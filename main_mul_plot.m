function main_mul_plot
% function main_mul_plot
%
%   Plots the results of main_mul_compute
%	from main_mul_results.mat
%    
%   Requires the ppfem package

%   Author: R. Andreev, Aug 2016
    
    suffix = '';
    
    D = load(['main_mul_results' suffix '.mat']);
    
    % Plotting nicer
    iif = @(C, varargin) varargin{2-logical(C)};
    nicer = iif(exist('./renice_plot.m', 'file'), @renice_plot, @(h)h);
    close all;
    
    fig_conv = figure;
    
    set(0, 'defaulttextinterpreter', 'latex');
    
    % Handle / method for legend
    HH = NaN(1, length(D.methods));
    for method = D.methods
        KK = D.KE{method == D.methods}.KK;
        EE = D.KE{method == D.methods}.EE;
    
        label = D.labels{method == D.methods};
        
        % Rate
        if (length(EE) >= 2)
            r = polyder(polyfit(log(KK), log(EE), 1));
            disp([label ' convergence rate: ' num2str(r)]);
        end
        disp('--');
        
        label = strrep(label, 'Q ', '');
        label = strrep(label, 'P ', '');
        label = strrep(label, 'iE', '{\textrm{iE}}');
        label = strrep(label, 'CN', '{\textrm{CN}}');
        label = strrep(label, '*', '{\star}');
        label = strrep(label, '\delta', 'box');
        label = strrep(label, 'box', '{\textrm{box}}');
        label = ['$' label '$']; % Assume latex interpreter
        D.labels{method == D.methods} = label;
        
        % Plot
        figure(fig_conv);
        HH(method == D.methods) = nicer(loglog(KK, EE, D.styles{method == D.methods}));
        hold on;
        title('Diagonal L1 error');
        xlabel('Temporal mesh width'); ylabel('Diagonal L1 error');
        pause(0.1);
    end
    
%     % Reverse order on the legend
%     reverse = @(X) X(length(X) : -1 : 1); HH = reverse(HH); D.labels = reverse(D.labels);

    h = nicer(legend(HH, D.labels{:}, 'Location', 'NW'));
    set(h, 'Interpreter', 'latex');
    
    % Save figure
    if (isequal(input('Save figure (y/n)? ', 's'), 'y'))
        figure(fig_conv);
        filename = ['conv_trace_error_M2' suffix];
        title(''); % No title in production figures
        saveas(gcf, filename, 'epsc');
        saveas(gcf, filename, 'png');
    end
end
