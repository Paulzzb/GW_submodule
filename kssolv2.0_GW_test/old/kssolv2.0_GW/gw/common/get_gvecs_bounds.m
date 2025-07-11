function [wfn_box_min, wfn_box_max] = get_gvecs_bounds(G, wfn_box_min, wfn_box_max)
G_min = double(min(G));
G_max = double(max(G));
wfn_box_min = min(G_min, wfn_box_min);
wfn_box_max = max(G_max, wfn_box_max);
end