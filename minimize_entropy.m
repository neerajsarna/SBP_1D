% we minize the entropy given the moments to find the maxwellian
function f = minimize_entropy(Ax,Ay,mass_matrix,rho,ux,uy,theta,id_sys)

vx = diag(Ax);
vy = diag(Ay);

value_f0 = arrayfun(@(x,y) f0(x,y),vx,vy);

switch id_sys
    case 1
        lagrange = mass_matrix\[ux; uy; rho; theta];
        f = [He1(vx),He1(vy),He0(vx),(He2(vx)+He2(vy))/sqrt(2)]*lagrange.*value_f0;        
    case 2
        lagrange = theta /(sqrt(2) * mass_matrix(3,3));
        f = lagrange * value_f0;
end


end

function f = He1(x)
f = x;
end

function f = He0(x)
f = ones(length(x),1);
end

function f = He2(x)
f = (x.^2 - 1)/sqrt(2);
end

% f0 corresponding to the 2d velocity space
function f = f0(vx,vy)
f = exp(-(vx^2+vy^2)/2)/(2 * pi);
end

