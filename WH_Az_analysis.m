%
% PURPOSE: Generate big ol' (not that big) matrix with 3 columns: (1) The q
% value (x in the plot) (2) the average standard deviation of peak width
% (3) the standard error in this value
%
% INPUTs: Row, Column, q values vector
%
%row = 27;
%col = 141;
function miso_Lcsd_err = WH_Az_analysis(row, col)

q_vect_short = [12 13 18 23 27];%[11 13 17 23 26];% 17
q_val = [0.01145 0.01308 0.01741 0.02272 0.02644];%0.0173 

num_fits = [5 6 6 5 6]; %5

q_w_err_mat = [];  %q_w_err is the initial output of peak-fitting script
 % This is Williamson Hall Plot values q^2 (um^-2), w^2 (um^-2), error (um^-2)

        % Finds plotting parameters and fits
        for i=1:length(q_vect_short)
           
            %datafile = sprintf('ZA-111/Col136_r25r65/All_Phi_LineProfs/%d_master00%d_ELp_0.0%d.dat',row,col,q_vect_short(i));
            datafile = sprintf('2018Jan_WHps/Dec17_WH_text/Azimuth_LineProfs_andCircInteg/_%d_master000%d_ELp_0.0%d.dat',row,col, q_vect_short(i));
            
            file = importdata(datafile,'\t',307);%,155);
            q_w_err_mat = [q_w_err_mat; AZ_width_fit(q_val(i), file.data, num_fits(i))];
        end
        q_w_err_mat;
%
%   IMPORTANT: IF Gaussian, exponent = 2 , if  Lorentzian, exponent = 1;
expo = 1.0;
%Calc the % distortion in the lattice spacing from the q_width_9avg vlaues
    q_locat = q_w_err_mat(:,1);  % inverse nm, vector of 5 q values
    q_width = q_w_err_mat(:,2);  % average width- also in inv nm
    wid_err = q_w_err_mat(:,3);  % error in width- also inv nm

    x_values = q_locat.^(expo); % in inv nm or inv square nanometer
    y_values = q_width.^(expo);
    y_errors = y_values.*expo.*wid_err./q_width; %y_values and q_width cancel
                                % out if expo = 1... wid_err doesn't change
    [param p_err] = linfitxy(x_values, y_values, zeros(length(x_values),1),y_errors) ;
    
   %{
   figure;
   if row<16  %36 before- just pick where color changes
       colr = 'green'; end
   if 15<row && row<31
       colr = 'red'; end
   if row>30
       colr = 'blue'; end
    errorbar([0;x_values], [param(2);y_values], [p_err(2);y_errors], 'o', 'LineWidth', 1.5,'MarkerSize',5,'Color',colr,'MarkerEdgeColor',colr,'MarkerFaceColor',colr);
    hold on
    plot(0:0.01:0.28, (0:0.01:0.28)*param(1)+param(2), colr);
    axis([-0.01 0.3 0 0.04]);
    xlabel('q (nm^{-1})'); ylabel('w_{\phi} (nm^{-1})');
    %title('Williamson Hall Plot of Azimuthal Peak Width'); 
     fig = gcf;
    pubgraph(fig, 24,1,'w')
    %}
    
    dist_or_miso = param(1)^(1/expo)*180/pi;
    d_or_m_err   = norm(dist_or_miso*(1/expo)*p_err(1)/param(1));
    Lcsd        =  ((2*pi)^(expo)/param(2))^(1/expo);
    Lcsd_err_neg  =  Lcsd - ((2*pi)^(expo)/(param(2)+p_err(2)))^(1/expo);%Lcsd*(1/expo)*(param(2)-p_err(2))/param(2)); 
    Lcsd_err_pos = ((2*pi)^(expo)/(param(2)-p_err(2)))^(1/expo) - Lcsd;

miso_Lcsd_err  = [dist_or_miso d_or_m_err Lcsd  Lcsd_err_neg Lcsd_err_pos];
