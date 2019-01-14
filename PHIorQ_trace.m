function [ steps_g_err_L_err ] = PHIorQ_trace( rows, stepsize, az_or_q )
%
% Rows vector should be the specific order ex: 45:-1:1 for data used in this study 
%   rows = 45:-1:1; 

% Stepsize is how far (in mm) each 2D SAXS image was acquired from the
% next. This is simply to set the x-axis of the trace. In this case:
%   stepsize = 0.1; % 0.1mm or 100 microns

% az_or_q  = 'az' % or
% az_or_q = 'q' % depending on whether the peak characteristics (specifically FWHM) 
% are measured along the azimuthal ('az') or scattering ('q') directions
%

col = 10;  % This is "down the center" of the map shown Fig S4c

stepvect = 0: stepsize : (length(rows)-1)*stepsize;
steps_g_err_L_err = zeros(length(rows),6);

if az_or_q == 'az'  %azimuthal analysis from cicular line profiles (5pxl width)
    Lcsd_lims = [0 7];
    miso_dist_lims = [0  10];
    for r = 1:length(rows)
     row_info = WH_Az_analysis(rows(r),col)
     steps_g_err_L_err(r,:) = [stepvect(r), row_info(1), row_info(2), row_info(3), row_info(4), row_info(5)]; 
    end

end

if az_or_q == 'q'
    Lcsd_lims = [0 40];
    miso_dist_lims = [0  0.1];
    for r = 1:length(rows)
        row_info = WH_q_analysis(rows(r),col)
        steps_g_err_L_err(r,:) = [stepvect(r), row_info(1), row_info(2), row_info(3), row_info(4), row_info(5)]; 
    end
end

 steps_g_err_L_err

% User must select whether to plot parameters based on WH graph slope
% (misorientation parameter or lattice distortion parameter) OR
% whether to plot data related to WH graph intercept
% (coherently scattering domain sizes)
%
%
%
        % Plot of g_phi (misorientation) or g_lat (Lattice Distortion (%))
%
figure
if az_or_q == 'az'
    errorbar(steps_g_err_L_err(:,1),steps_g_err_L_err(:,2),steps_g_err_L_err(:,3),'s','Color','black')
    axis([-0.1 4.6 miso_dist_lims(1) miso_dist_lims(2)]);
    ylabel('g_{\phi} (deg)'); 
end
if az_or_q == 'q'  % Lattice distortion data must be converted from fraction to %
    errorbar(steps_g_err_L_err(:,1),steps_g_err_L_err(:,2)*100,steps_g_err_L_err(:,3)*100,'s','Color','black')
    axis([-0.1 4.6 miso_dist_lims(1)*100 miso_dist_lims(2)*100]);
   ylabel('Lattice Distortion (%)')
end
xlabel('Line scan location (mm)');
fig = gcf;
pubgraph(fig, 24,1,'white')
%}
 
%{
    % Plot of L_csd (
figure
errorbar(steps_g_err_L_err(:,1),steps_g_err_L_err(:,4)/10^3,steps_g_err_L_err(:,5)/10^3,steps_g_err_L_err(:,6)/10^3,'s','Color','black')
axis([-0.1 4.6 Lcsd_lims(1) Lcsd_lims(2)]);
    xlabel('Line scan location (mm)'); ylabel('L_{CSD} (\mum)');
    %title('Williamson Hall Plot of q Peak Width'); 
     fig = gcf;
    pubgraph(fig, 24,1,'white')
    %}
