% we develop the 
function B_ID1 = dvlp_B_ID1(varargin)

B_ID2 = varargin{1};

if nargin == 1
    % number of columns are the number of variables
    id_odd = 2:2:size(B_ID2,2);
end

if nargin == 2
    [id_odd,~] = get_id_Odd(varargin{2});
    id_odd = flatten_cell(id_odd);
end

    B_ID1 = B_ID2;
    
    % flip the signs of the odd variables
    B_ID1(:,id_odd) = -B_ID1(:,id_odd);
    
end