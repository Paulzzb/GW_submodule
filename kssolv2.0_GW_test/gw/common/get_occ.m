function [occ] = get_occ(options, ik, ispin)
ev = options.ev(:, ik, ispin);
efermi = options.efermi;
tol = 0.001;
for i = 1:size(ev,1)
    if (ev(i) > (efermi + tol))
        occ(i) = 0;
    elseif (ev(i) < (efermi - tol))
        occ(i) = 1;
    else
        occ(i) = 0.5;
    end
end