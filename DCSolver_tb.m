OPT.LUT_PATH = 'D:\New folder\TSMC_180_05';
if exist('NCHV_ac','var') == 0
        load([OPT.LUT_PATH, '/NCHV_ac.mat'])
end

if exist('PCHV_ac','var') == 0
    load([OPT.LUT_PATH, '/PCHV_ac.mat'])
end
models.nmos = NCHV_ac;
models.pmos = PCHV_ac;
maximum_step = 50;
[Circuit, x] = DCSolver('netlist_5TOTA_1.txt', models, maximum_step);