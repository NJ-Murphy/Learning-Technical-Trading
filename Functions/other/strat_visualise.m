function outputs = strat_visualise(varargin)
% strat_visualise calls the w-th strategy function and retruns the buy/sell signals
% for the w-th strategy at time t. The function will also return any other
% information required by each strategy such as the value of the strategy
% at the previous time increment. 
%
% strat_visualise(xt, strategyfn, p.ntype,ell1,k1)
%
% strat_visualise(xt, strategyfn, p.ntype,ell1)

%% initialise the input variables
xt = varargin{1};
strategyfn = varargin{2};
Inputs = varargin{3}';
ell1 = varargin{4};
optargs = {xt,strategyfn,Inputs,ell1,{}};
% now put these defaults into the valuesToUse cell array
optargs(1:nargin) = varargin;

clear varargin

% Place optional args in memorable variable names 
[xt,strategyfn,Inputs,ell1,k1] = optargs{:};


%% Compute the wealth of the w-th expert/strategy

% Find what inputs the strategy function requires
for i = 1:length(Inputs)
   if strfind(Inputs{i},'n2')~=0 
       stratargs{1} = xt;
       stratargs{2} = ell1;
   else
       stratargs{1} = xt;
       stratargs{2} = ell1;
       stratargs{3} = k1; 
   end
end

% Call the w-th strategy and compute the signals
outputs = strategyfn(stratargs(:));

end