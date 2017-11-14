function printFig(varargin)
% printFig  Print a figure window as eps file. 
%   printFig(path)  write current figure to path
%   printFig(figHandle,path)   figHandle is a figure or axes handle

if nargin == 2
    fig = varargin{1};
    if isa(fig,'matlab.graphics.axis.Axes')
        fig = fig.Parent;
    end
    path = varargin{2};
else
    fig = gcf;
    path = varargin{1};
end
print(fig,path,'-depsc');
end