function Circuit_desc = parser(netlist)
% Parser function supports only MOSFETs, Resistors, indep. sources
% Parser is case-insenstive
% All nodes should be non-negative integers 1,2,3...
% Ground node name: 0
% MOSFETS: Mxx pmos/nmos model_name drain_node gate_node source_node bulk_node
%            [ VG0=value VD0=value VS0=value VB0=value L=value W=value ID=value RHO=value ]
% Resistors: Rxx node+ node- [ r=value i=value v=value ]
% Indep. voltage sources: Vxx node+ node- [v=value i0=value]
% Indep current sources: Ixx node+ node- [i=value]
fid = fopen(netlist,'rt');
if fid < 0
    error('error opening netlist.txt');
end

mosfets_i = 1;
linear_i = 1;
max_node = 0;
no_v_sources = 0;
no_of_mosfets_to_design = 0;
while true
    oneline = lower(fgetl(fid));
    if ~ischar(oneline)
        break;
    end
    line_splitted = strsplit(strtrim(oneline));
    identifier = line_splitted{1}(1);
    if identifier == 'm'
        Circuit_desc.mosfets{mosfets_i}.name = line_splitted{1};
        Circuit_desc.mosfets{mosfets_i}.type = line_splitted{2};
        Circuit_desc.mosfets{mosfets_i}.model = line_splitted{3};
        if ~strcmp(line_splitted{3},'solve') 
            no_of_mosfets_to_design = no_of_mosfets_to_design + 1;
        end
        Circuit_desc.mosfets{mosfets_i}.drain = str2double(line_splitted{4});
        Circuit_desc.mosfets{mosfets_i}.gate = str2double(line_splitted{5});
        Circuit_desc.mosfets{mosfets_i}.source = str2double(line_splitted{6});
        Circuit_desc.mosfets{mosfets_i}.bulk = str2double(line_splitted{7});
        for j=8:length(line_splitted)
            par_value = strsplit(line_splitted{j},'=');
            Circuit_desc.mosfets{mosfets_i}.(par_value{1}) = str2double(par_value{2});
        end
        max_node_i = max(Circuit_desc.mosfets{mosfets_i}.drain,...
            max(Circuit_desc.mosfets{mosfets_i}.gate,...
            max(Circuit_desc.mosfets{mosfets_i}.source,...
            Circuit_desc.mosfets{mosfets_i}.bulk)));
        if max_node_i > max_node
            max_node = max_node_i;
        end
        mosfets_i = mosfets_i + 1;
    else
        Circuit_desc.linear{linear_i}.name = line_splitted{1};
        Circuit_desc.linear{linear_i}.pnode = str2double(line_splitted{2});
        Circuit_desc.linear{linear_i}.nnode = str2double(line_splitted{3});
        for j=4:length(line_splitted)
            par_value = strsplit(line_splitted{j},'=');
            Circuit_desc.linear{linear_i}.(par_value{1}) = str2double(par_value{2});
        end
        max_node_i = max(Circuit_desc.linear{linear_i}.pnode,...
            Circuit_desc.linear{linear_i}.nnode);
        if max_node_i > max_node
            max_node = max_node_i;
        end
        if identifier == 'v'
            no_v_sources = no_v_sources + 1;
        end
        linear_i = linear_i + 1;
    end
end
Circuit_desc.no_of_nodes = max_node;
Circuit_desc.no_of_vsources = no_v_sources;
Circuit_desc.no_of_design_mosfets = no_of_mosfets_to_design;