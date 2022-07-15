function [G_linearized, RHS_linearized, H, g] = mosfet_stamper(Circuit_desc)
dim = Circuit_desc.no_of_nodes + Circuit_desc.no_of_vsources ...
    + Circuit_desc.no_of_design_mosfets;
no_of_mosfets = length(Circuit_desc.mosfets);
H = zeros(dim,no_of_mosfets);
G_linearized = zeros(dim,dim);
RHS_linearized = zeros(dim,1);
i = 1;
g = zeros(no_of_mosfets,1);
i_design_mosfets = Circuit_desc.no_of_nodes + Circuit_desc.no_of_vsources + 1;
for element = Circuit_desc.mosfets
    if element{1}.drain ~= 0
        H(element{1}.drain,i) = 1;
    end
    if element{1}.source ~= 0
        H(element{1}.source,i) = -1;
    end
    if strcmp(element{1}.model,'solve')
        if element{1}.drain ~= 0
            G_linearized(element{1}.drain,element{1}.drain) = ...
                G_linearized(element{1}.drain,element{1}.drain) + element{1}.gds;
        end

        if element{1}.drain ~= 0 && element{1}.source ~= 0
            G_linearized(element{1}.source,element{1}.drain) = ...
                G_linearized(element{1}.source,element{1}.drain) - element{1}.gds;

            G_linearized(element{1}.drain,element{1}.source) = ...
                G_linearized(element{1}.drain,element{1}.source) ...
                - element{1}.gds - element{1}.gm + element{1}.gmb;
        end

        if element{1}.source ~= 0
            G_linearized(element{1}.source,element{1}.source) = ...
                G_linearized(element{1}.source,element{1}.source) ...
                + element{1}.gds + element{1}.gm - element{1}.gmb;
        end

        if element{1}.drain ~= 0 && element{1}.gate ~= 0
            G_linearized(element{1}.drain,element{1}.gate) = ...
                G_linearized(element{1}.drain,element{1}.gate) + element{1}.gm;
        end

        if element{1}.gate ~= 0 && element{1}.source ~= 0
            G_linearized(element{1}.source,element{1}.gate) = ...
                G_linearized(element{1}.source,element{1}.gate) - element{1}.gm;
        end

        if element{1}.gate ~= 0
            G_linearized(element{1}.gate,element{1}.gate) = ...
                G_linearized(element{1}.gate,element{1}.gate);
        end

        if element{1}.drain ~= 0 && element{1}.bulk ~= 0
            G_linearized(element{1}.drain,element{1}.bulk) = ...
                G_linearized(element{1}.drain,element{1}.bulk) - element{1}.gmb;
        end

        if element{1}.source ~= 0 && element{1}.bulk ~= 0
            G_linearized(element{1}.source,element{1}.bulk) = ...
                G_linearized(element{1}.source,element{1}.bulk) + element{1}.gmb;
        end

        if element{1}.drain ~= 0
            RHS_linearized(element{1}.drain) = ...
                RHS_linearized(element{1}.drain) -  element{1}.Ieq;
        end

        if element{1}.source ~= 0
            RHS_linearized(element{1}.source) = ...
                RHS_linearized(element{1}.source) +  element{1}.Ieq;
        end
        g(i) = element{1}.ID;
    elseif strcmp(element{1}.model,'design')
        if element{1}.drain ~= 0
            G_linearized(i_design_mosfets, element{1}.drain) = ...
                G_linearized(i_design_mosfets,element{1}.drain) - element{1}.gds/element{1}.gm;
            RHS_linearized(element{1}.drain) = ...
                RHS_linearized(element{1}.drain) -  element{1}.ID;
        end

        if element{1}.source ~= 0
            G_linearized(i_design_mosfets, element{1}.source) = ...
                G_linearized(i_design_mosfets,element{1}.source) - 1 + (element{1}.gds - element{1}.gmb)/element{1}.gm;
            
            G_linearized(element{1}.source,i_design_mosfets) = ...
                G_linearized(element{1}.source, i_design_mosfets) - 1;
            RHS_linearized(element{1}.source) = ...
                RHS_linearized(element{1}.source) +  element{1}.ID;
        end
        
        if element{1}.gate ~= 0
            G_linearized(element{1}.gate,i_design_mosfets) = ...
                G_linearized(element{1}.gate, i_design_mosfets) + 1;
            G_linearized(i_design_mosfets, element{1}.gate) = ...
                G_linearized(i_design_mosfets,element{1}.gate) + 1;
            RHS_linearized(i_design_mosfets) = ...
                RHS_linearized(i_design_mosfets) +  element{1}.VGSeq;
        end
        
        if element{1}.bulk ~= 0
            G_linearized(i_design_mosfets, element{1}.bulk) = ...
                G_linearized(i_design_mosfets,element{1}.bulk) + element{1}.gmb/element{1}.gm;
        end
        i_design_mosfets = i_design_mosfets + 1;
        g(i) = element{1}.id;
    end
    g(i) = element{1}.ID;
    i = i + 1;
end