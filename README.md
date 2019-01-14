# Williamson_Hall_Analysis

Data hierarchy:

% In order to create a graph showing the spatial variation of a Williamson-Hall parameter (a "trace"):
PHIorQ_trace( rows, stepsize, az_or_q ) % must be run. Sample settings for processing the data in this study are:
%
% Rows vector should be the specific order ex: 45:-1:1 for data used in this study 
%   rows = 45:-1:1; Each element points towards a set of datafiles extracted from the 2D SAXS patterns using manuscript reference 47

% Stepsize is how far (in mm) each 2D SAXS image was acquired from the
% next. This is simply to set the x-axis of the trace. In this case:
%   stepsize = 0.1; % 0.1mm or 100 microns

% az_or_q  = 'az' % or
% az_or_q = 'q' % depending on whether the peak characteristics (specifically FWHM) 
% are measured along the azimuthal ('az') or scattering ('q') directions
%

PHIorQ_trace() will iterate through all rows. At each row it has the dependencies on:

WH_Az_analysis(rows(r),col)  % if 'az' is set for az_or_q
WH_q_analysis(rows(r),col)  % if 'q' is set for az_or_q

Each of the two scripts WH_Az_analysis and WH_q_analysis will compile all the peak fitting data for a given row, and column
As such, they depend on the peak fitting algorithms to fit each individual Bragg peak:

AZ_width_fit() or 
Q_width_fit() depending on if the peaks are being fit according to their azimuthal (az) or scattering (q) 
  direction profiles.
  
  The raw text data extracted from the 2D SAXS using profile fitting from manuscipt reference 47 is given in the compressed data file.
  It should 2018Jan_WHps.zip should be extracted prior to running PHIorQ_trace.
