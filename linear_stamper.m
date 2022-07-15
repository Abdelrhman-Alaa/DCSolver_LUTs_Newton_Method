function [G, RHS, x0] = linear_stamper(Circuit_desc, x0)
dim = Circuit_desc.no_of_nodes + Circuit_desc.no_of_vsources ...
    + Circuit_desc.no_of_design_mosfets;
G = zeros(dim,dim);
RHS = zeros(dim,1);
vsrc_i = 1;
for element = Circuit_desc.linear
    if element{1}.name(1) == 'r'
        if element{1}.pnode ~= 0
            G(element{1}.pnode,element{1}.pnode) = ...
                G(element{1}.pnode,element{1}.pnode) + 1/element{1}.r;
        end
        if element{1}.nnode ~= 0
            G(element{1}.nnode,element{1}.nnode) = ...
                G(element{1}.nnode,element{1}.nnode) + 1/element{1}.r;
        end
        if element{1}.nnode ~= 0 && element{1}.nnode ~= 0
            G(element{1}.pnode,element{1}.nnode) = ...
                G(element{1}.pnode,element{1}.nnode) - 1/element{1}.r;
            G(element{1}.nnode,element{1}.pnode) = ...
                G(element{1}.nnode,element{1}.pnode) - 1/element{1}.r;
        end
    end
    if element{1}.name(1) == 'i'
        if element{1}.pnode ~= 0
            RHS(element{1}.pnode) = ...
                RHS(element{1}.pnode) - element{1}.i;
        end
        if element{1}.nnode ~= 0
            RHS(element{1}.nnode) = ...
                RHS(element{1}.nnode) + element{1}.i;
        end
    end
    if element{1}.name(1) == 'v'
        if element{1}.pnode ~= 0
            G(element{1}.pnode,Circuit_desc.no_of_nodes + vsrc_i) = 1;
            G(Circuit_desc.no_of_nodes + vsrc_i,element{1}.pnode) = 1;
        end
        if element{1}.nnode ~= 0
            G(element{1}.nnode,Circuit_desc.no_of_nodes + vsrc_i) = -1;
            G(Circuit_desc.no_of_nodes + vsrc_i,element{1}.nnode) = -1;
        else
            x0(element{1}.pnode) = element{1}.v;
        end
        RHS(Circuit_desc.no_of_nodes + vsrc_i) = element{1}.v;
        vsrc_i = vsrc_i + 1;
    end
end