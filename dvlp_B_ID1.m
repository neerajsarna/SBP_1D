% we develop the 
function B_ID1 = dvlp_B_ID1(B_ID2)

% number of columns are the number of variables
id_odd = 2:2:size(B_ID2,2);

B_ID1 = B_ID2;

% flip the signs of the odd variables
B_ID1(:,id_odd) = -B_ID1(:,id_odd);
end