%% Bolus Vax1 Changing C0
p = base_case_parameters();
[txt_name, writeoption] = getName('Vax1_01_C0scan');
C0s = p.C0 * [8, 4, 2, 1, 1/2, 1/4, 1/8];
for C0 = C0s
    p.C0 = C0;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax1 Scan - Varying E1h and dE12
p = base_case_parameters();
[txt_name, writeoption] = getName('Vax1_02_scan_E1h_dE12');
for E1h = 7:0.2:8
    for dE12 = 0:0.2:1
        p.E1h = E1h; p.dE12 = dE12;
        createParameters(txt_name, writeoption, p)
    end
end
%% Bolus Vax1 Scan - Varying dE12 and p
p = base_case_parameters();
[txt_name, writeoption] = getName('Vax1_03_scan_dE12_and_p');
for dE12 = [0, 0.2, 0.4, 0.6, 0.8, 1]
    for subp = [0.5,0.25,0.1,0.02]
        p.dE12 = dE12;
        p.p = subp;
        createParameters(txt_name, writeoption, p)
    end
end
%% Bolus Vax1 No Epitope Masking
p = base_case_parameters();
[txt_name, writeoption] = getName('Vax1_04_nomasking');
createParameters(txt_name, writeoption, p)
%% Bolus Vax1 Epitope Masking, with/without steric hindrance
p = base_case_parameters();
p.masking = 1;
[txt_name, writeoption] = getName('Vax1_05_masking');
createParameters(txt_name, writeoption, p)
%Changing steric hindrance
[txt_name, writeoption] = getName('Vax1_06_steric_change');
for steric = [0.2, 0.3, 0.4, 0.5, 0.7]
    p.steric = steric;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax2 Changing C0
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 180;
[txt_name, writeoption] = getName('Vax2_01_C0scan');
C0s = p.C0 * [8, 4, 2, 1, 1/2, 1/4, 1/8];
for C0 = C0s
    p.C0 = C0;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax2 Scan - Varying E1h and dE12
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 180;
[txt_name, writeoption] = getName('Vax2_02_scan_E1h_dE12');
for E1h = 7:0.2:8
    for dE12 = 0:0.2:1
        p.E1h = E1h; p.dE12 = dE12;
        createParameters(txt_name, writeoption, p)
    end
end
%% Bolus Vax2 Scan - Varying dE12 and p
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 180;
[txt_name, writeoption] = getName('Vax2_03_scan_dE12_and_p');
for dE12 = [0, 0.2, 0.4, 0.6, 0.8, 1]
    for subp = [0.5,0.25,0.1,0.02]
        p.dE12 = dE12;
        p.p = subp;
        createParameters(txt_name, writeoption, p)
    end
end
%% Bolus Vax2 No Epitope Masking
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 180;
[txt_name, writeoption] = getName('Vax2_04_nomasking');
createParameters(txt_name, writeoption, p)

%% Bolus Vax2 Epitope Masking, with/without steric hindrance
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 180;
p.masking = 1;
[txt_name, writeoption] = getName('Vax2_05_masking');
createParameters(txt_name, writeoption, p)
%Changing steric hindrance
[txt_name, writeoption] = getName('Vax2_06_steric_change');
for steric = [0.2, 0.3, 0.4, 0.5, 0.7]
    p.steric = steric;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax2 42 days
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 42;
[txt_name, writeoption] = getName('Vax2_07_42d');
createParameters(txt_name, writeoption, p)
p.masking = 1;
p.steric = 0.3;
createParameters(txt_name, writeoption, p)
%% Bolus Vax2 changing memToGCFraction
p = base_case_parameters();
p.vaxnum = 2;
p.tmax = 180;
[txt_name, writeoption] = getName('Vax2_08_memToGCFraction');
for memToGC = [0.01, 0.02, 0.04, 0.08]
    p.memToGCFrac = memToGC;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax3 No Epitope Masking
p = base_case_parameters();
p.vaxnum = 3;
p.tmax = 28;
[txt_name, writeoption] = getName('Vax3_01_nomasking');
createParameters(txt_name, writeoption, p)

%% Bolus Vax3 Epitope Masking, with/without steric hindrance
p = base_case_parameters();
p.vaxnum = 3;
p.tmax = 28;
p.masking = 1;
[txt_name, writeoption] = getName('Vax3_02_masking');
createParameters(txt_name, writeoption, p)
%Changing steric hindrance
[txt_name, writeoption] = getName('Vax3_03_steric_0.3');
for steric = [0.3]
    p.steric = steric;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax3 Long, with/without Masking Steric Hindrance
p = base_case_parameters();
p.vaxnum = 3;
p.tmax = 180;
[txt_name, writeoption] = getName('Vax3_04_180d');
createParameters(txt_name, writeoption, p)
%Changing steric hindrance
p.masking = 1;
for steric = [0.3]
    p.steric = steric;
    createParameters(txt_name, writeoption, p)
end
%% Bolus Vax3 Early Booster
p = base_case_parameters();
p.vaxnum = 3;
p.tmax = 28;
p.earlybooster = 42;
[txt_name, writeoption] = getName('Vax3_05_earlybooster');
createParameters(txt_name, writeoption, p)
p.masking = 1;
p.steric = 0.3;
createParameters(txt_name, writeoption, p)

%% Bolus Vax4, Masking Steric Hindrance
p = base_case_parameters();
p.vaxnum = 4;
p.tmax = 28;
[txt_name, writeoption] = getName('Vax4');
createParameters(txt_name, writeoption, p)
%Changing steric hindrance
p.masking = 1;
for steric = [0.3]
    p.steric = steric;
    createParameters(txt_name, writeoption, p)
end
%% Changing rho
p = base_case_parameters();
[txt_name, writeoption] = getName('ChangingRho');
for rho=[0.9, 0.97]
    p.rho = rho;
    createParameters(txt_name, writeoption, p, 1)
end
%% Changing output 
p = base_case_parameters();
[txt_name, writeoption] = getName('ChangingOutput');
for outputprob = [0.05, 0.10]
    for outputpcfrac = [0.1, 0.5]
    p.outputprob = outputprob;
    p.outputpcfrac = outputpcfrac;
    createParameters(txt_name, writeoption, p, 1)
    end
end
%% Changing w2
p = base_case_parameters();
[txt_name, writeoption] = getName('Changingw2');
for w2=[0.3, 0.7, 1]
    p.w2 = w2;
    createParameters(txt_name, writeoption, p, 1)
end
%% Changing w1
p = base_case_parameters();
p.w2 = 1;
[txt_name, writeoption] = getName('Changingw1');
for w1=[20, 100, 500]
    p.w1 = w1;
    createParameters(txt_name, writeoption, p, 1)
end

%% Changing C0
p = base_case_parameters();
[txt_name, writeoption] = getName('Changing_C0');
C0s = p.C0 * [8, 4, 2, 1, 1/2, 1/4, 1/8];
for C0 = C0s
    p.C0 = C0;
    createParameters(txt_name, writeoption, p, 1)
end

%% Bolus Vax1 Scan - Varying E1h and dE12
p = base_case_parameters();
[txt_name, writeoption] = getName('Changing_E1h_dE12');
for E1h = 7:0.2:8
    for dE12 = 0:0.2:1
        p.E1h = E1h; p.dE12 = dE12;
        createParameters(txt_name, writeoption, p, 1)
    end
end
%% Bolus Vax1 Scan - Varying dE12 and p
p = base_case_parameters();
[txt_name, writeoption] = getName('Changing_dE12_and_p');
for dE12 = [0, 0.2, 0.4, 0.6, 0.8, 1]
    for subp = [0.5,0.25,0.1,0.02]
        p.dE12 = dE12;
        p.p = subp;
        createParameters(txt_name, writeoption, p, 1)
    end
end
