function out=makepositive(in)
if isa(in,'int16')
    % Typecast causes negative int16 values to overflow, adding 2^16
    % Typecast only works on flat arrays, so first flatten the data
    positive_values=typecast(in(:), 'uint16');

    % Convert data to single and unflatten
    out=reshape(single(positive_values), size(in));
else
    out=in;
end
end