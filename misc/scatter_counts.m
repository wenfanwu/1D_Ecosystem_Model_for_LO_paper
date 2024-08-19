function [counts, xbins, ybins] = scatter_counts(obs_val, sim_val, bin_width)
% Prepare for scatter-density diagram
%
%% Syntax
% [counts, xbins, ybins] = scatter_counts(obs_val, sim_val)
% [counts, xbins, ybins] = scatter_counts(obs_val, sim_val, bin_width)
%
%% Description
% [counts, xbins, ybins] = scatter_counts(obs_val, sim_val) prepares the
% data for scatter-density diagrams.
% [counts, xbins, ybins] = scatter_counts(obs_val, sim_val, bin_width)
% specifies the bin width along x- and y- dimensions.
%
%% Examples 
% clc;clearvars
% obs_val = randn(1,1000000);
% sim_val = obs_val + randn(1,1000000)*0.5;
% [counts, xbins, ybins] = scatter_counts(obs_val, sim_val)
% figure; pcolor(xbins, ybins, counts);shading interp
%
%% Input Arguments
% obs_val --- observed data
% sim_val --- simulated data
% bin_width --- width of bins along each dimension (x and y). e.g. bin_width = [0.1 0.1];
%
%% Output Arguments
% counts --- bin counts, specified as a matrix. 
% xbins --- bin sequence in x-dimension.
% ybins --- bin sequence in y-dimension.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institue of Marine Science in 2024. 
% Last Updated on 2024-8-18. 
% Email: wwu@vims.edu
% 
% See also: histcounts

%% Parse inputs
if nargin==2
    obs_gaps = (max(obs_val, [], 'omitnan')-min(obs_val, [], 'omitnan'))/100;
    sim_gaps = (max(sim_val, [], 'omitnan')-min(sim_val, [], 'omitnan'))/100;
    min_gap =  min(obs_gaps,sim_gaps);
    bin_width = [min_gap min_gap];
end

%% Calculate
[counts, xbins, ybins] = histcounts2(obs_val, sim_val,'BinWidth',bin_width);

xbins = xbins + bin_width(1);
ybins = ybins + bin_width(2);
xbins(end) = [];
ybins(end) = [];

counts(counts==0) = nan;
counts = counts';

end