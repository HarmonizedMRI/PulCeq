function s = setfields(s, varargin)
% setfields - set one or more struct fields
%
% setfields(s,'field1',value1,'field2',value2,...)
%
% Replace one or more fields in a struct

for k = 1:2:numel(varargin)
    s.(varargin{k}) = varargin{k+1};
end

