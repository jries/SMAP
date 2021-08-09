function copysize(fig, ref, tar, posRef, posTar)
    % fig: figure obj
    % ref: the tag of the reference
    % tar: the tag of the target
    % posRef and posTar: position index
    h_ref = findobj(fig, 'tag', ref);
    h_tar = findobj(fig, 'tag', tar);
    set([h_ref h_tar], 'Units', 'Centimeter')
    h_tar.Position(posTar) = h_ref.Position(posRef);
    set([h_ref h_tar], 'Units', 'Normalized')
end