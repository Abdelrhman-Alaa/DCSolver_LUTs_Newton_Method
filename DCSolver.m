function [Circuit, x] = DCSolver(netlist, models, maximum_step, relative_tol, vabsolute_tol, iabsolute_tol)
if nargin > 2
    max_step = maximum_step;
else
    max_step = 10;
end
if nargin > 3
    reltol = relative_tol;
else
    reltol = 0.01;
end
if nargin > 4
    vabstol = vabsolute_tol;
else
    vabstol = 1e-6;
end
if nargin > 5
    iabstol = iabsolute_tol;
else
    iabstol = 1e-9;
end
Circuit = parser(netlist);
[Circuit, x0] = dev_eval(Circuit,models);
[G, RHS, x0] = linear_stamper(Circuit, x0);
x = cell(1,20);
x_linear = cell(1,20);
x{1} = x0;
ord_of_conv = zeros(1,20);
for i = 2:20
    [G_linearized, RHS_linearized, H, g] = mosfet_stamper(Circuit);
    Js = G + G_linearized;
    RHS_total = RHS + RHS_linearized;
    x_linear{i} = Js \ RHS_total;
    bool_big_step = abs(x_linear{i} - x{i - 1}) > max_step;
    if i == 2
        bool_big_step = bool_big_step & (x{1} > eps);
        Residual = G*x{1} + H*g - RHS;
        Residual_i_1 = norm(Residual(1:Circuit.no_of_nodes));
        Residual_v_1 = norm(Residual(Circuit.no_of_nodes+1:Circuit.no_of_nodes+Circuit.no_of_vsources));
    end
    x{i} = x_linear{i} .* (1 - bool_big_step) + ...
        bool_big_step .* (x{i - 1} + sign(x_linear{i} - x{i - 1}) .* max_step);

    dx_v = norm(x{i}(1:Circuit.no_of_nodes) - x{i-1}(1:Circuit.no_of_nodes));
    dx_i = norm(x{i}(Circuit.no_of_nodes+1:end) ...
        - x{i-1}(Circuit.no_of_nodes+1:end));
    if norm(x{i}) > norm(x{i-1})
        xmax = x{i};
    else
        xmax = x{i-1};
    end
    xmax_v = norm(xmax(1:Circuit.no_of_nodes));
    xmax_i = norm(xmax(Circuit.no_of_nodes+1:Circuit.no_of_nodes+Circuit.no_of_vsources));
    Residual = G*x{i} + H*g - RHS;
    Residual_i = norm(Residual(1:Circuit.no_of_nodes));
    Residual_v = norm(Residual(Circuit.no_of_nodes+1:Circuit.no_of_nodes+Circuit.no_of_vsources));
    if i > 3
        ord_of_conv(i - 3) = log(norm(x{i} - x{i-1})/norm(x{i-1} - x{i-2})) ...
            / log(norm(x{i-1} - x{i-2})/norm(x{i-2} - x{i-3}));
    end
    if dx_v < reltol * xmax_v + vabstol && dx_i < reltol * xmax_i + iabstol ...
            && Residual_v < reltol * Residual_v_1 + vabstol ...
            && Residual_i < reltol * Residual_i_1 + iabstol 
       Circuit = dev_eval(Circuit,models,x{i}); 
       break;
    end
    Circuit = dev_eval(Circuit,models,x{i});
end