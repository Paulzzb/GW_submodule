function [gcutoff] = gcutoff(ng, ekin, isrtrq, ecutoff)
if (max(ekin) < ecutoff)
    gcutoff = ng;
else
    gup = ng;
    gdn = 1;
    for ig = 1:ng
        gmid = round((gup + gdn)/2);
        if (gmid == gdn)
            break
        end
        if (ekin(isrtrq(gmid)) > ecutoff)
            gup = gmid;
        else
            gdn = gmid;
        end
    end
    gcutoff = gdn;
end
