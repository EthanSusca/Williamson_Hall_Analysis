%
%
%  Will loop through many files at every q value in order to accumulated
%  statistics on the width in the scattering vector direction as well as
%  the more accurate value of the scattering vector to begin with
%
function LattD_Lcsd_LattDerr_Lcsderr = WH_q_analysis(row,col)

%bunch of vectors corresponding to equivalent slices
%27_141   0.24 % +/-0.03.
q211type = [58 118 178 238 298 ] ; %129 missing(249) 69:60:309; % 5 (211) and (422) points
q220type = [28 88 148 208 268 328];%219 %missing 159 %39:60:339; % 6 (220) and (440) points
q321type = [47 69 167 189 227 249 ];% [everythign else]... all these locations +/- 8 degrees] 
%57_137 =  0.24 % +/0.02 for 40_139 = .18+/-.02
%q211type = [9 69 129 189  ];%309] ; % missing(249) 69:60:309; % 5 (211) and (422) points
%q220type = [39 99 159 219 279];% 339];%missing 159 %39:60:339; % 6 (220) and (440) points
%error_adjust


%q_vect = [0.011392  0.013129 0.0173  0.0229 0.026425];% 0.0174
q_vect = [0.0115 0.0134 0.0176 0.0231 0.0264];

vect211 =[]; %zeros(length(q211type),2);
vect422 =[]; %zeros(length(q211type),2);
vect321 = [];
vect220 =[]; %zeros(length(q220type),2);
vect440 =[]; %zeros(length(q220type),2);
%Nested for loops for evaluating each q value, of which there are two in
%each file
%FIRST DEALL WITH 211 and 422, then 220 and 440
   % figure
   % hold on
    
for q = 1:length(q_vect)
    if round(q_vect(q),3)==0.012%0.011  % this is {211}
        for id = 1:length(q211type)
            datafile = sprintf('2018Jan_WHps/Dec17_WH_text/Qwidth_sector/_%d_master000%d_%d_8.dat',row,col,q211type(id));
            file = importdata(datafile,'\t',304);
            vect211 = [vect211; Q_width_fit(file.data, q_vect(q))];
        end
    end
    if round(q_vect(q),3)==0.013  % this is {220}
        for id = 1:length(q220type)
            datafile = sprintf('2018Jan_WHps/Dec17_WH_text/Qwidth_sector/_%d_master000%d_%d_8.dat',row,col,q220type(id));
            file = importdata(datafile,'\t',304);
            vect220 = [vect220; Q_width_fit(file.data, q_vect(q))];
        end
    end
    if round(q_vect(q),3)==0.018  % this is {321}
        for id = 1:length(q321type)
            datafile = sprintf('2018Jan_WHps/Dec17_WH_text/Qwidth_sector/_%d_master000%d_%d_6.dat',row,col,q321type(id));
            file = importdata(datafile,'\t',307);
            vect321 = [vect321; Q_width_fit(file.data, q_vect(q))];
        end
    end
    if round(q_vect(q),3)==0.023  % this is {422}
        for id = 1:length(q211type)
            datafile = sprintf('2018Jan_WHps/Dec17_WH_text/Qwidth_sector/_%d_master000%d_%d_8.dat',row,col,q211type(id));
            file = importdata(datafile,'\t',304);
            vect422 = [vect422; Q_width_fit(file.data, q_vect(q))];
        end
    end
    if round(q_vect(q),3)==0.026  % this is {440}
        for id = 1:length(q220type)
            datafile = sprintf('2018Jan_WHps/Dec17_WH_text/Qwidth_sector/_%d_master000%d_%d_8.dat',row,col,q220type(id));
            file = importdata(datafile,'\t',304);
            vect440 = [vect440; Q_width_fit(file.data, q_vect(q))];
        end
    end
end
%hold off

vect211;
vect422;
vect220;
vect440;
vect321;

q211 = [mean(vect211(:,1)) mean(vect211(:,2)) std(vect211(:,2))/(length(vect211(:,2))^0.5)];
q220 = [mean(vect220(:,1)) mean(vect220(:,2)) std(vect220(:,2))/(length(vect220(:,2))^0.5)];
q321 = [mean(vect321(:,1)) mean(vect321(:,2)) std(vect321(:,2))/(length(vect321(:,2))^0.5)];
q422 = [mean(vect422(:,1)) mean(vect422(:,2)) std(vect422(:,2))/(length(vect422(:,2))^0.5)];
q440 = [mean(vect440(:,1)) mean(vect440(:,2)) std(vect440(:,2))/(length(vect440(:,2))^0.5)];

preWHmat = [q211; q220; q321; q422; q440];%q321;
    expo = 1;
    %Calc the % distortion in the lattice spacing from the q_width_9avg vlaues
    q_locat = preWHmat(:,1);  % inverse nm, vector of 5 q values
    q_width = preWHmat(:,2);  % average width- also in inv nm
    wid_err = preWHmat(:,3);  % error in width- also inv nm

    WH_mat(:,1) = q_locat.^(expo); % in inv nm or inv square nanometer
    WH_mat(:,2) = q_width.^(expo);
    WH_mat(:,3) = WH_mat(:,2).*expo.*wid_err./q_width; %y_values and q_width cancel
                                % out if expo = 1... wid_err doesn't change
%     
%     [param p_err] = linfitxy(x_values, y_values, zeros(length(WH_mat),1),WH_mat(:,3)) ;

[param, p_err] = linfitxy(WH_mat(:,1),WH_mat(:,2),zeros(length(WH_mat(:,3)),1),WH_mat(:,3));
%hold off
 %{
 %figure;
    if row<16  %36 before- just pick where color changes
       colr = 'green'; end
   if 15<row && row<31
       colr = 'red'; end
   if row>30
       colr = 'blue'; end
    errorbar([0;WH_mat(:,1)], [param(2);WH_mat(:,2)], [p_err(2);WH_mat(:,3)],'o', 'LineWidth', 1.5,'MarkerSize',5,'Color',colr,'MarkerEdgeColor',colr,'MarkerFaceColor',colr);%, 's', 'LineWidth', 1.5);
    hold on
    plot(0:0.01:0.28, (0:0.01:0.28)*param(1)+param(2), 'b');
    axis([-0.01 0.3 -0.005 0.015]);
    xlabel('q_{q} (nm^{-1})'); ylabel('w_{q} (nm^{-1})');
    title('WH Plots: FWHM along scattering vector (q)'); 
    fig = gcf;
    pubgraph(fig, 24,1,'white')
    %}
       %
   %figure;

    
 

LatDist = (param(1))^(1/expo); % in degrees
LatDist_err = norm(LatDist*(1/expo)*p_err(1)/param(1));
% COHERENT SCATTERING DOMAIN LENGTH and ERROR
Lcsd = ((2*pi)^(expo)/param(2))^(1/expo);
%Lcsd_err = norm(Lcsd*(1/expo)*p_err(2)/param(2));
Lcsd_err_neg  =  Lcsd - ((2*pi)^(expo)/(param(2)+p_err(2)))^(1/expo);%Lcsd*(1/expo)*(param(2)-p_err(2))/param(2)); 
Lcsd_err_pos = ((2*pi)^(expo)/(param(2)-p_err(2)))^(1/expo) - Lcsd;

% RETURN 1x4 vector
format shortG;
LattD_Lcsd_LattDerr_Lcsderr  = [LatDist  LatDist_err Lcsd Lcsd_err_neg Lcsd_err_pos];



%%%
%%  OLD, for the record:

            
%         
% %     
% % for n = 1:length(q211type)
% %     for m = 1:length(q_vect)
%         datafile = sprintf('ZA-111/try1_10degsect/%d_master00%d_%d_5.dat',row,col,q211type(n));
%         file = importdata(datafile,'\t',150);
%         if m == 1
%             vect211 = [vect211; Q_width_fit(file.data, q_vect(m))];%, qloc211(m)*0.1, 10^7)];
%         end
%         if m == 2
%             vect422 = [vect422; Q_width_fit(file.data, q_vect(m))];%, qloc211(m)*0.1, 5*10^5)];
%         end
%     end
% end
%     hold off
% % Now the same thing, but for the 220/440 set
%     figure
%     hold on
% for n = 1:length(q220type)
%     for m = 1:length(replacemeqloc220)
%         datafile = sprintf('ZA-111/try1_10degsect/%d_master00%d_%d_5.dat',row,col,q220type(n));
%         file = importdata(datafile,'\t',150);
%         if m == 1
%             vect220 = [vect220; Q_width_fit(file.data, replacemeqloc220(m))];%, qloc220(m)*0.1, 10^6)];
%         end
%         if m == 2
%             vect440 = [vect440; Q_width_fit(file.data, replacemeqloc220(m))];%, qloc220(m)*0.1, 8*10^4)];
%         end
%     end
% end
% hold off

        
        
        
        
        
        