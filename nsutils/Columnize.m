function columnVector = Columnize(tensor)
    columnVector = reshape(tensor,[numel(tensor) 1]);
end