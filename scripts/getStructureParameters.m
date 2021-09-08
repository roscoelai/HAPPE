function p = getStructureParameters(mystruct, myfield, default_value)
% Returns mystruct.myfield if it exists, or default_value otherwise
if isfield(mystruct, myfield)
    p = mystruct.(myfield);
elseif exist("default_value", "var")
    p = default_value;
else
    error("Either 'mystruct.myfield' or 'default_value' must exist");
end
