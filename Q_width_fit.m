function [ q_invnm_width ] = Q_width_fit( data, q_guess)%, width_guess, offset_guess )
            %  INPUTS
% datafile
% a string that is the address of the image to analyze

% q_guess 
%  guess the mean of the distribution we are going to fit-
% width_guess
%   guess a HWHM value- fit more often successful if larger than actual HWHM

            %  OUTPUTS
% q_invnm       % the location of the average across all successful fits

% width     % the actual HWHM... the gamma term in a Lorentz
% wid_err   % the standard error in the widths found

% n_count    % just return a 1 if the fit is successful (q determined to be
                % a reasoanble value +/- 0.001 of q_invnm_guess

%file = importdata(datafile,'\t',36);

% Look for max over a certain range defined by +/- some # data points 

search_range = 5;%5;  %+/- to look for max point
[val qind] = min(abs(data(:,1)-q_guess));
[mx ndx] =  max(data((qind-search_range):(qind+search_range),2));
indx = find(data(:,2) == mx);

fit_range = 9; % how many points +/- around max to make fit to
X = data(indx-fit_range:indx+fit_range,1);
Y = data(indx-fit_range:indx+fit_range,2);
data(indx(1),1);

p0 = [0.1, data(indx(1),1), data(indx(1),1)*0.015,  0];%mx*0.0005];%0.01,5];

% p0 = [multiplier(scalar)  x-mean-mu-value  width-about-x  y-offset]
%LORENTZIAN FIT
lorentz = @(p) norm( p(1)./(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2))+p(4)-Y(:));
[p resid] = fminsearchbnd(lorentz, p0, [0 p0(2)-0.001 0.0001 0], [mx*10 p0(2)+0.001 0.001 mx]);%...
    %*(10^6* (p(1) < 0)|(p(3)<0)|(p(4) < 0)), p0);

% GAUSSIAN FIT
%[p resid] = fminsearch(@(p) norm( p(1).*exp(-0.5.*((X(:)-p(2))./p(3)).^2 )./((2*pi).^0.5.*p(3))+p(4) -Y(:)), p0);
% PSUEDO VOIGT
%p0 = [p0, 0.5]; %add another fit term that is 0 for Lorentz and 1 for Gauss Distributions
%[p resid] = fminsearch(@(p) norm(  p(1) .*...
  %  (   p(5).*exp(-0.5* ((X(:)-p(2))./(2.*p(3)/2.355)).^2)./((2*pi).^0.5.*p(3)) +...
   %   (1-p(5)).*(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2)).^(-1) )+  p(4)  -  Y(:) ), p0);

%p0
%p
%plot(X, log(Y))
%plot(X,log(p(1)./(pi.*p(3).*(1+((X(:)-p(2))./p(3)).^2))+p(4)))


if abs(p(2)-q_guess)<0.0005 && (p(3) < 0.002)
    q_invnm_width = [10*p(2), 10*norm(2*p(3))];%, p(5)] ;

   
else
    q_invnm_width = [] ;
end

             
                


end

