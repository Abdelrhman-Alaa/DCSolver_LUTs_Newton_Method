function [Circuit_desc, x0] = dev_eval(Circuit_desc,models,x)
i = 1;
dim = Circuit_desc.no_of_nodes + Circuit_desc.no_of_vsources ...
    + Circuit_desc.no_of_design_mosfets;
x0 = zeros(dim, 1);
for element = Circuit_desc.mosfets
    if nargin == 2
        if isfield(element{1}, 'vg0')
            Circuit_desc.mosfets{i}.VG = element{1}.vg0;
            if element{1}.gate ~= 0
                x0(element{1}.gate) = element{1}.vg0;
            end
        else
            Circuit_desc.mosfets{i}.VG = 0;
        end
        if isfield(element{1}, 'vd0')
            Circuit_desc.mosfets{i}.VD = element{1}.vd0;
            if element{1}.drain ~= 0
                x0(element{1}.drain) = element{1}.vd0;
            end
        else
            Circuit_desc.mosfets{i}.VD = 0;
        end
        if isfield(element{1}, 'vs0')
            Circuit_desc.mosfets{i}.VS = element{1}.vs0;
            if element{1}.source ~= 0
                x0(element{1}.source) = element{1}.vs0;
            end
        else
            Circuit_desc.mosfets{i}.VS = 0;
        end
        if isfield(element{1}, 'vb0')
            Circuit_desc.mosfets{i}.VB = element{1}.vb0;
            if element{1}.bulk ~= 0
                x0(element{1}.bulk) = element{1}.vb0;
            end
        else
            Circuit_desc.mosfets{i}.VB = 0;
        end
    else
        if element{1}.gate ~= 0
            Circuit_desc.mosfets{i}.VG = x(element{1}.gate);
        end
        if element{1}.drain ~= 0
            Circuit_desc.mosfets{i}.VD = x(element{1}.drain);
        end
        if element{1}.source ~= 0
            Circuit_desc.mosfets{i}.VS = x(element{1}.source);
        end
        if element{1}.bulk ~= 0
            Circuit_desc.mosfets{i}.VB = x(element{1}.bulk);
        end
    end
    if strcmp(element{1}.type,'nmos')
        type_flag = 1;
    elseif strcmp(element{1}.type,'pmos')
        type_flag = -1;
    end
    if strcmp(element{1}.model,'solve')
        L = element{1}.l * 1e6;
        if isfield(element{1},'varw')
            Circuit_desc.mosfets{i}.w = Circuit_desc.mosfets{Circuit_desc.mosfets{i}.varw}.W * 1e-6;
            if nargin == 2
                Circuit_desc.mosfets{i}.VGS = Circuit_desc.mosfets{Circuit_desc.mosfets{i}.varw}.VGS;
                VGS = Circuit_desc.mosfets{i}.VGS;
            else
                VGS = Circuit_desc.mosfets{i}.VG - Circuit_desc.mosfets{i}.VS;
            end
        else
            VGS = Circuit_desc.mosfets{i}.VG - Circuit_desc.mosfets{i}.VS;
        end
        W = Circuit_desc.mosfets{i}.w * 1e6;
        VDS = Circuit_desc.mosfets{i}.VD - Circuit_desc.mosfets{i}.VS;
        VSB = Circuit_desc.mosfets{i}.VS - Circuit_desc.mosfets{i}.VB;
        TEMP = 27;
                
        Circuit_desc.mosfets{i}.ID = type_flag*aaLookupID(models.(element{1}.type), L,...
            VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP) / models.(element{1}.type).W * W;
        Circuit_desc.mosfets{i}.gm = aaLookupGM(models.(element{1}.type), L,...
            VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP) / models.(element{1}.type).W * W;
        Circuit_desc.mosfets{i}.gds = ...
            models.(element{1}.type).GDS(L, VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP)...
            / models.(element{1}.type).W * W;
        Circuit_desc.mosfets{i}.gmb = models.(element{1}.type).GMB(L, VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP)...
            / models.(element{1}.type).W * W;
        Circuit_desc.mosfets{i}.Ieq = Circuit_desc.mosfets{i}.ID - ...
            Circuit_desc.mosfets{i}.gm * VGS - ...
            Circuit_desc.mosfets{i}.gds * VDS - ...
            Circuit_desc.mosfets{i}.gmb * VSB;
    else
        L = element{1}.l * 1e6;
        RHO = element{1}.rho;
        VDS = Circuit_desc.mosfets{i}.VD - Circuit_desc.mosfets{i}.VS;
        VSB = Circuit_desc.mosfets{i}.VS - Circuit_desc.mosfets{i}.VB;
        TEMP = 27;
        Circuit_desc.mosfets{i}.VGS = type_flag*models.(element{1}.type).VGS_RHO(L,...
            RHO, VDS*type_flag, VSB*type_flag, TEMP);
        VGS = Circuit_desc.mosfets{i}.VGS;
        Circuit_desc.mosfets{i}.VG = Circuit_desc.mosfets{i}.VS ...
            + Circuit_desc.mosfets{i}.VGS;
        % This part should be enhanced. The value of W is not needed to
        % be calculated each iteration. normalized values for gm, gds, gmb
        % should be used. The value of W should be calculated at last
        % iteration and then compute the final values of gm, gds, gmb
        if nargin == 2
            W = Circuit_desc.mosfets{i}.id / aaLookupID(models.(element{1}.type), L,...
            VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP) * models.(element{1}.type).W;
        else
            W = element{1}.W;
        end
        Circuit_desc.mosfets{i}.ID0_W = aaLookupID(models.(element{1}.type), L,...
            VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP) / models.(element{1}.type).W;
        Circuit_desc.mosfets{i}.gm = aaLookupGM(models.(element{1}.type), L,...
            VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP) / models.(element{1}.type).W * W;
        Circuit_desc.mosfets{i}.gds = ...
            models.(element{1}.type).GDS(L, VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP)...
            / models.(element{1}.type).W * W;
        Circuit_desc.mosfets{i}.gmb = models.(element{1}.type).GMB(L, VGS*type_flag, VDS*type_flag, VSB*type_flag, TEMP)...
            / models.(element{1}.type).W * W;
        
        Circuit_desc.mosfets{i}.VGSeq = Circuit_desc.mosfets{i}.VGS - ...
            Circuit_desc.mosfets{i}.gds / Circuit_desc.mosfets{i}.gm * VDS - ...
            Circuit_desc.mosfets{i}.gmb / Circuit_desc.mosfets{i}.gm * VSB;
        
        Circuit_desc.mosfets{i}.W = Circuit_desc.mosfets{i}.id ...
            / Circuit_desc.mosfets{i}.ID0_W;
        
        Circuit_desc.mosfets{i}.ID = Circuit_desc.mosfets{i}.id*type_flag;
    end
    i = i + 1;
end