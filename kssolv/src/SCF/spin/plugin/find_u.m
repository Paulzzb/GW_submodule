function u = find_u(s)
 %  Receive as input a 3*3 rotation matrix s, and gives as output 
 % the matrix u which represents the same rotation in the spin space
 % only useful when using symmetry.
if abs(det(s)+1) < 1e-8
    saux = -s;
else
    saux = s;
end

if sum(abs(saux-eye(3)) >= 1e-8,'all') == 0
    u = eye(2);
end
end