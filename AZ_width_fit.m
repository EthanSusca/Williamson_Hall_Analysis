function [ q_w_err ] = AZ_width_fit( qval_i_AtoNm, data, n_fits )
% qval_invnm (must be in inverse nanomaters
% data: fist column degrees, second is intensity
% n_fits: how many distributions to fit with single azimuthal line profile

q_w_err = [];
w_invnm_frm_deg = [];
for n = 1: n_fits
        [mx, indx] =  max(data(:,5));
        data;  mx;  indx;
       
        %make indx-30 or +30 more reasonable value
        rang = indx-24:indx+24; %so +/-6 degrees
        
        if any(rang(:)<1) || any(rang(:)>length(data))
            p = [1 2 3] ; %forcing p(3) < 30 means don't add datapoint to list
        else
            X = data(rang,4); % data is fit to +/- 32/4 = 8 degrees about a peak
            Y = data(rang,5);
            p0 = [mx, data(indx(1),4), 5, 0.01*mx ];%mx*0.005];%0.01,5];
            % p0 = [multiplier(scalar)  x-mean-mu-value  width-about-x  y-offset]
            %LORENTZIAN FIT
            [p resid] = fminsearch(@(p) norm( p(1)./(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2))+p(4)-Y(:)), p0);
            data(rang,5) = data(rang,5)- (p(1)./(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2))+p(4));
            %plot(X,log10(p(1)./(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2))+p(4)));
            % GAUSSIAN FIT
            %[p resid] = fminsearch(@(p) norm( p(1).*exp(-0.5.*((X(:)-p(2))./p(3)).^2 )./((2*pi).^0.5.*p(3))+p(4) -Y(:)), p0);
            %data(rang,5) = data(rang,5)- (p(1).*exp(-0.5.*((X(:)-p(2))./p(3)).^2 )./((2*pi).^0.5.*p(3))+p(4));
            % PSUEDO VOIGT
            %p0 = [p0, 0.5]; %add another fit term that is 0 for Lorentz and 1 for Gauss Distributions
            %[p resid] = fminsearch(@(p) norm(  p(1) .*...
              %  (   p(5).*exp(-0.5* ((X(:)-p(2))./(2.*p(3)/2.355)).^2)./((2*pi).^0.5.*p(3)) +...
               %   (1-p(5)).*(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2)).^(-1) )+  p(4)  -  Y(:) ), p0);
            %plot(X,log10(p(1)./(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2)))
        end     
     
     if p(3) < 30
        w_invnm_frm_deg = [w_invnm_frm_deg; 2*p(3)*pi/180*qval_i_AtoNm*10]; %arc length in 1/nm       
     else
         w_invnm_frm_deg = [w_invnm_frm_deg;[]];
     end
end

    q_w_err = [qval_i_AtoNm*10, mean(w_invnm_frm_deg), std(w_invnm_frm_deg)/(n_fits)^0.5, n_fits];%, p(5)] ;



end

