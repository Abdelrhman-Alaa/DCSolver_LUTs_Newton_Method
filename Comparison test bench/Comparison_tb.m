%M1 nmos design 2 1 0 0 L=1 VD0=0.2 id=10 rho=15
%M2 nmos design 4 3 2 0 L=1 VD0=0.4 VS0=0.2 VB0=0 rho=15 id=10
%M3 pmos design 4 4 5 5 L=1 VD0=0.6 VS0=1 VB0=1 rho=15 id=10
%v2 3 0 v=0.7
%v4 5 0 v=1
OPT.LUT_PATH = 'D:\New folder\TSMC_180_05';
if exist('NCHV_ac','var') == 0
        load([OPT.LUT_PATH, '/NCHV_ac.mat'])
end
if exist('PCHV_ac','var') == 0
    load([OPT.LUT_PATH, '/PCHV_ac.mat'])
end
L1_vec = [1 3 4 6 8];
RHO1_vec = [10 12 14 16 18];
L2_vec = [1 3 4 6 8];
RHO2_vec = [10 12 14 16 18];
L3_vec = [1 3 4 6 8];
RHO3_vec = [10 12 14 16 18];
TEMP = 27;
v2 = 0.7;
v4 = 1;
t_newton = zeros(1,6);
newton_flags = ones(1,6);
t_fixed_point = zeros(1,6);
fixed_point_flags = ones(1,6);
tic
for L1 = L1_vec
    for L2 = L2_vec
        for L3 = L3_vec
            for RHO1 = RHO1_vec
                for RHO2 = RHO2_vec
                    for RHO3 = RHO3_vec
                        %L1 = 1;
                        VDS1 = 0.2;
                        VSB1 = 0;
                        %RHO1 = 15;
                        ID1 = 10e-6;
                        VGS1 = NCHV_ac.VGS_RHO(L1, RHO1, VDS1, VSB1, TEMP);
                        gm1 = aaLookupGM(NCHV_ac, L1, VGS1, VDS1, VSB1, TEMP);
                        gds1 = NCHV_ac.GDS(L1, VGS1, VDS1, VSB1, TEMP);
                        gmb1 = NCHV_ac.GMB(L1, VGS1, VDS1, VSB1, TEMP);
                        VGSeq1 = VGS1 - gds1 / gm1 * VDS1 - gmb1 / gm1 * VSB1;

                        %L2 = 1;
                        VDS2 = 0.2;
                        VSB2 = 0.2;
                        %RHO2 = 15;
                        ID2 = 10e-6;
                        VGS2 = NCHV_ac.VGS_RHO(L2, RHO2, VDS2, VSB2, TEMP);
                        gm2 = aaLookupGM(NCHV_ac, L2, VGS2, VDS2, VSB2, TEMP);
                        gds2 = NCHV_ac.GDS(L2, VGS2, VDS2, VSB2, TEMP);
                        gmb2 = NCHV_ac.GMB(L2, VGS2, VDS2, VSB2, TEMP);
                        VGSeq2 = VGS2 - gds2 / gm2 * VDS2 - gmb2 / gm2 * VSB2;

                        %L3 = 1;
                        VDS3 = 0.4;
                        VSB3 = 0;
                        %RHO3 = 15;
                        ID3 = 10e-6;
                        TEMP = 27;
                        VGS3 = PCHV_ac.VGS_RHO(L3, RHO3, VDS3, VSB3, TEMP);
                        gm3 = aaLookupGM(PCHV_ac, L3, VGS3, VDS3, VSB3, TEMP);
                        gds3 = PCHV_ac.GDS(L3, VGS3, VDS3, VSB3, TEMP);
                        gmb3 = PCHV_ac.GMB(L3, VGS3, VDS3, VSB3, TEMP);
                        VGSeq3 = (VGS3 - gds3 / gm3 * VDS3 - gmb3 / gm3 * VSB3) * -1;
                        x = cell(1,20);
                        for i=1:20  
                              Js = [ 0,                     0, 0,            0,                                0, 0, 0, 1,  0,  0 ...
                                ; 0,                     0, 0,            0,                                0, 0, 0, 0, -1,  0 ...
                                ; 0,                     0, 0,            0,                                0, 1, 0, 0,  1,  0 ...
                                ; 0,                     0, 0,            0,                                0, 0, 0, 0,  0,  1 ...
                                ; 0,                     0, 0,            0,                                0, 0, 1, 0,  0, -1 ...
                                ; 0,                     0, 1,            0,                                0, 0, 0, 0,  0,  0 ...
                                ; 0,                     0, 0,            0,                                1, 0, 0, 0,  0,  0 ...
                                ; 1,             -gds1/gm1, 0,            0,                                0, 0, 0, 0,  0,  0 ...
                                ; 0, (gds2 - gmb2)/gm2 - 1, 1,    -gds2/gm2,                                0, 0, 0, 0,  0,  0 ...
                                ; 0,                     0, 0, 1 - gds3/gm3, gmb3/gm3 + (gds3 - gmb3)/gm3 - 1, 0, 0, 0,  0,  0];
                              RHS_total = [...
                                       0; ...
                                       0; ...
                                       0; ...
                                       0; ...
                                     -ID3; ...
                                      v2; ...
                                      v4; ...
                                  VGSeq1; ...
                                  VGSeq2; ...
                                  VGSeq3];
                                x{i} = Js \ RHS_total;
                                VDS1 = x{i}(2);
                                VDS2 = x{i}(4) - x{i}(2);
                                VSB2 = x{i}(2);
                                VDS3 = x{i}(5) - x{i}(4);

                                VGS1 = NCHV_ac.VGS_RHO(L1, RHO1, VDS1, VSB1, TEMP);
                                gm1 = aaLookupGM(NCHV_ac, L1, VGS1, VDS1, VSB1, TEMP);
                                gds1 = NCHV_ac.GDS(L1, VGS1, VDS1, VSB1, TEMP);
                                gmb1 = NCHV_ac.GMB(L1, VGS1, VDS1, VSB1, TEMP);
                                VGSeq1 = VGS1 - gds1 / gm1 * VDS1 - gmb1 / gm1 * VSB1;

                                VGS2 = NCHV_ac.VGS_RHO(L2, RHO2, VDS2, VSB2, TEMP);
                                gm2 = aaLookupGM(NCHV_ac, L2, VGS2, VDS2, VSB2, TEMP);
                                gds2 = NCHV_ac.GDS(L2, VGS2, VDS2, VSB2, TEMP);
                                gmb2 = NCHV_ac.GMB(L2, VGS2, VDS2, VSB2, TEMP);
                                VGSeq2 = VGS2 - gds2 / gm2 * VDS2 - gmb2 / gm2 * VSB2;

                                VGS3 = PCHV_ac.VGS_RHO(L3, RHO3, VDS3, VSB3, TEMP);
                                gm3 = aaLookupGM(PCHV_ac, L3, VGS3, VDS3, VSB3, TEMP);
                                gds3 = PCHV_ac.GDS(L3, VGS3, VDS3, VSB3, TEMP);
                                gmb3 = PCHV_ac.GMB(L3, VGS3, VDS3, VSB3, TEMP);
                                VGSeq3 = (VGS3 - gds3 / gm3 * VDS3 - gmb3 / gm3 * VSB3) * -1;
                                
                                if i > 1
                                    dx_v = norm(x{i}(1:5) - x{i-1}(1:5));
                                    if dx_v < 0.01 * norm(x{i - 1}(1:5))
                                        break
                                    end
                                end
                        end
                    end
                    if newton_flags(1) == 1
                        t_newton(1) = toc;
                        newton_flags(1) = 0;
                    end
                end
                if newton_flags(2) == 1
                    t_newton(2) = toc;
                    newton_flags(2) = 0;
                end
            end
            if newton_flags(3) == 1
                t_newton(3) = toc;
                newton_flags(3) = 0;
            end
        end
        if newton_flags(4) == 1
            t_newton(4) = toc;
            newton_flags(4) = 0;
        end
    end
    if newton_flags(5) == 1
        t_newton(5) = toc;
        newton_flags(5) = 0;
    end
end
t_newton(6) = toc;
tic
ID0_W1 = aaLookupID(NCHV_ac, L1, VGS1, VDS1, VSB1, TEMP) / NCHV_ac.W;
W1 = ID1 / ID0_W1;
gm1 = gm1 / NCHV_ac.W * W1;
gds1 = gds1 / NCHV_ac.W * W1;
gmb1 = gmb1 / NCHV_ac.W * W1;

ID0_W2 = aaLookupID(NCHV_ac, L2, VGS2, VDS2, VSB2, TEMP) / NCHV_ac.W;
W2 = ID2 / ID0_W2;
gm2 = gm2 / NCHV_ac.W * W2;
gds2 = gds2 / NCHV_ac.W * W2;
gmb2 = gmb2 / NCHV_ac.W * W2;

ID0_W3 = aaLookupID(PCHV_ac, L3, VGS3, VDS3, VSB3, TEMP) / NCHV_ac.W;
W3 = ID3 / ID0_W3;
gm3 = gm3 / NCHV_ac.W * W3;
gds3 = gds3 / NCHV_ac.W * W3;
gmb3 = gmb3 / NCHV_ac.W * W3;
t_to_add = toc;
t_newton = t_newton + t_to_add;
tic
for L1 = L1_vec
    for L2 = L2_vec
        for L3 = L3_vec
            for RHO1 = RHO1_vec
                for RHO2 = RHO2_vec
                    for RHO3 = RHO3_vec
                        VDS3 = 0.4;
                        VSB3 = 0;
                        for k = 1:20
                            VGS3 = PCHV_ac.VGS_RHO(L3, RHO3, VDS3, VSB3, TEMP);
                            if abs(VGS3 - VDS3) / VDS3 < 0.01 * VDS3
                                break
                            end
                            VDS3 = VGS3;
                        end
                        VDS3 = VGS3;
                        VDS2 = v4 - VGS3 - 0.2;
                        VSB2 = 0.2;
                        for k = 1:20
                            VGS2 = NCHV_ac.VGS_RHO(L2, RHO2, VDS2, VSB2, TEMP);
                            VDS2_new = v4 - VGS3 - 0.7 + VGS2;
                            VSB2_new = 0.7 - VGS2;
                            if abs(VDS2_new - VDS2) / VDS2 < 0.01 * VDS2 ...
                                    && abs(VSB2_new - VSB2) / VSB2 < 0.01 * VSB2                                    
                                break
                            end
                            VDS2 = VDS2_new;
                            VSB2 = VSB2_new;
                        end
                        VDS2 = VDS2_new;
                        VSB2 = VSB2_new;
                        VDS1 = VSB2;
                        VSB1 = 0;
                        VGS1 = NCHV_ac.VGS_RHO(L1, RHO1, VDS1, VSB1, TEMP);
                    end
                    if fixed_point_flags(1) == 1
                        t_fixed_point(1) = toc;
                        fixed_point_flags(1) = 0;
                    end
                end
                if fixed_point_flags(2) == 1
                    t_fixed_point(2) = toc;
                    fixed_point_flags(2) = 0;
                end
            end
            if fixed_point_flags(3) == 1
                t_fixed_point(3) = toc;
                fixed_point_flags(3) = 0;
            end
        end
        if fixed_point_flags(4) == 1
            t_fixed_point(4) = toc;
            fixed_point_flags(4) = 0;
        end
    end
    if fixed_point_flags(5) == 1
        t_fixed_point(5) = toc;
        fixed_point_flags(5) = 0;
    end
end
t_fixed_point(6) = toc;
tic
ID0_W1 = aaLookupID(NCHV_ac, L1, VGS1, VDS1, VSB1, TEMP) / NCHV_ac.W;
W1 = ID1 / ID0_W1;
gm1 = gm1 / NCHV_ac.W * W1;
gds1 = gds1 / NCHV_ac.W * W1;
gmb1 = gmb1 / NCHV_ac.W * W1;

ID0_W2 = aaLookupID(NCHV_ac, L2, VGS2, VDS2, VSB2, TEMP) / NCHV_ac.W;
W2 = ID2 / ID0_W2;
gm2 = gm2 / NCHV_ac.W * W2;
gds2 = gds2 / NCHV_ac.W * W2;
gmb2 = gmb2 / NCHV_ac.W * W2;

ID0_W3 = aaLookupID(PCHV_ac, L3, VGS3, VDS3, VSB3, TEMP) / NCHV_ac.W;
W3 = ID3 / ID0_W3;
gm3 = gm3 / NCHV_ac.W * W3;
gds3 = gds3 / NCHV_ac.W * W3;
gmb3 = gmb3 / NCHV_ac.W * W3;
t_to_add = toc;
for i = 1:6
    t_fixed_point(i) = t_fixed_point(i) + t_to_add;
end
number_of_circuits = 6.^(1:6);
semilogx(number_of_circuits, t_fixed_point, number_of_circuits, t_newton);
semilogx(number_of_circuits, t_newton ./ t_fixed_point);
figure(1);
semilogx(number_of_circuits, t_fixed_point, number_of_circuits, t_newton);
legend('Fixed-Point time','Newton Method time');
xlabel('Number of Design Points');
ylabel('Time (Seconds)');
grid on;
figure(2);
semilogx(number_of_circuits, t_newton ./ t_fixed_point);
xlabel('Number of Design Points');
ylabel('Newton Method time / Fixed-Point time');
grid on;