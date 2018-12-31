function var = fluxVarGrabber(F,varpath,catdim)
% function var = fluxVarGrabber(F,varpath,catdim)
% Drills down into a multi-leg flight flux structure and extracts a specified variable for all legs.
%
% INPUTS:
% F:        flux structure. Should contain sub-structures for each leg ('L1','L2',etc).
% varpath:  cell array specifying names of sub-fields and, lastly, the variable of interest.
% catdim:   OPTIONAL scalar specifying dimension to concatenate variables from each leg.
%           Default is to use the longest dimension of the variable.
%
% OUTPUTS:
% var:      concatenated matrix of requested variables from each of the legs.
%
% EXAMPLE USE:
% varpath = {'T','wave','flux'};
% Tflux = fluxVarGrabber(F,varpath,1);
%
% 20161004 GMW

%get legs
Legs = fieldnames(F);
isaleg = strncmp('L',Legs,1);
Legs(~isaleg) = [];
nLegs = length(Legs);

% get variables
var=cell(nLegs,1);
for i=1:nLegs
    Lvar = F.(Legs{i});
    for j=1:length(varpath)
        if ~isfield(Lvar,varpath{j})
            disp(['fluxVarGrabber: field "' varpath{j} ': not found.'])
            var = [];
            return
        else
            Lvar = Lvar.(varpath{j}); %drill down
        end
    end
    var{i} = Lvar;
end

% determine cat dimension if not specified
if nargin<3 || isempty(catdim)
    s = size(var{1});
    [~,catdim] = max(s); %use longest dimension
end

% concatenate
var = cat(catdim,var{:});


