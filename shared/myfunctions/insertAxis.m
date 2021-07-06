function insertAxis(f, axTarget, tag_axRef)
    axRef = findobj(f, 'tag',tag_axRef);
    axTarget.Parent = f;
    axTarget.Units = 'Normalized';
    axRef.Units = 'Normalized';
    axTarget.Position = axRef.Position;
    axTarget.Tag = axRef.Tag;
    delete(axRef);
end