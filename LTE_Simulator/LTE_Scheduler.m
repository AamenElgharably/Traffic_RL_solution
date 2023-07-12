clear all;
close all;
clc;

%% Global Variables

Ga_Tx = 10;                 % Tx antenna gain (dB)
Ga_Rx = 5;                  % Rx antenna gain (dB)
Num_UEs = 50;               % Number of UEs

NoiseFigure = randi([5 10]);            % Its range = [5, 10] dB       
NoisePower_dBm = -174 + 10*log10(20*10^6) + NoiseFigure;
% -174 dBm is the power spectral density of thermal noise at room temperature (T = 293 K)

%% Configuration of the three eNBs
% All eNBs are aligned in the x direction, spaced by 35000m.
% All distances are in meters.

Spacing = 35e3;
for m = 1:3
    eNB(m).x = 5e3 + Spacing*(m-1);
    eNB(m).y = 25e3;
    eNB(m).Height = 30;
    eNB(m).Bandwidth = 20;                  % 20 MHz
    eNB(m).CarrierFreq = 1500;              % 1500 MHz
    eNB(m).PowerTransmit = 30;              % dBm
    eNB(m).ResourceBlocks = 100;            % Number of RBs in one slot
end

%% Allocation of the UEs

for i = 1:Num_UEs
    UE(i).x = randi([(eNB(2).x) - 25e3, (eNB(2).x) + 25e3]);
    UE(i).y = randi([(eNB(2).y) - 13e3, (eNB(2).y) + 13e3]);
    UE(i).z = 1.5;
end

for i = 1:Num_UEs
    Cord_x(i) = (UE(i).x)./1e3;
    Cord_y(i) = (UE(i).y)./1e3;
end

Cord_eNBx(1) = (eNB(1).x)./1e3;
Cord_eNBx(2) = (eNB(2).x)./1e3;
Cord_eNBx(3) = (eNB(3).x)./1e3;
Cord_eNBy(1:3) = (eNB(m).y)./1e3;

%% Path Loss and Power Receive Calculations

for m = 1:3
    for i = 1:Num_UEs
        % Path Loss Calculation
        UE(i).Distance(m) = Distance_Calculation(eNB(m).x, eNB(m).y, eNB(m).Height, UE(i).x, UE(i).y, UE(i).z);
        UE(i).PathLoss_dB(m) = (Okumura_HATA_Model((UE(i).Distance(m)./1e3), eNB(m).CarrierFreq, eNB(m).Height, UE(i).z));        
        % Power Receive Calculation
        UE(i).PowerReceive_dBm(m) = eNB(m).PowerTransmit - UE(i).PathLoss_dB(m) + Ga_Tx + Ga_Rx;
    end
end

%% Handover Parameters: RSSI, RSRP, RSRQ

% RSSI: Received Signal Strength Indicator
% RSRP: Reference Signal Received Power
% RSRQ: Reference Signal Received Quality
% Online Calculator 1: https://comtech.vsb.cz/qualmob/rsrq.html
% Online Calculator 2: https://arimas.com/2017/11/06/lte-rsrp-rsrq-rssi-calculator/

for i = 1:Num_UEs
    Pr1(i) = UE(i).PowerReceive_dBm(1);
    Pr2(i) = UE(i).PowerReceive_dBm(2);
    Pr3(i) = UE(i).PowerReceive_dBm(3);
    UE(i).RSSI_dBm(1) = Pr1(i) + NoisePower_dBm - (Pr2(i) + Pr3(i));
    UE(i).RSSI_dBm(2) = Pr2(i) + NoisePower_dBm - (Pr1(i) + Pr3(i));
    UE(i).RSSI_dBm(3) = Pr3(i) + NoisePower_dBm - (Pr1(i) + Pr2(i));
end

for m = 1:3    
    for i = 1:Num_UEs
        UE(i).RSRP_dBm(m) = UE(i).RSSI_dBm(m) - (10*log10(12*eNB(m).ResourceBlocks));
        % Range = [-140, -44] dBm
        UE(i).RSRQ_dB(m) = 10*log10(eNB(m).ResourceBlocks) + UE(i).RSRP_dBm(m) - UE(i).RSSI_dBm(m);
        % Range = [-19.5, -3] dB
    end
end

%% Handover Decision (Hysteresis-based Algorithm)

Hys = 3;      % Hysteresis offset, its Range = [0, 10] dB
eNB(1).CIO = 5; eNB(2).CIO = 0; eNB(3).CIO = 5;

for i = 1:Num_UEs
    % We assume that the initial serving eNB is the middle one (eNB 2)
    UE(i).ServingeNB = 2;
	% Check handover condition for eNB 1
	if UE(i).RSRP_dBm(1) + eNB(1).CIO > UE(i).RSRP_dBm(2) + eNB(2).CIO + Hys
        UE(i).ServingeNB = 1;
	% Check handover condition for eNB 3    
    elseif UE(i).RSRP_dBm(3) + eNB(3).CIO > UE(i).RSRP_dBm(2) + eNB(2).CIO + Hys
        UE(i).ServingeNB = 3;
    end	
end

%% SINR Calculation

for i = 1:Num_UEs
	if UE(i).ServingeNB == 1
        UE(i).Interference_dBm = UE(i).PowerReceive_dBm(2) + UE(i).PowerReceive_dBm(3);
        UE(i).PowerReceive_dBm = UE(i).PowerReceive_dBm(1);
    elseif UE(i).ServingeNB == 2
        UE(i).Interference_dBm = UE(i).PowerReceive_dBm(1) + UE(i).PowerReceive_dBm(3);
        UE(i).PowerReceive_dBm = UE(i).PowerReceive_dBm(2);
    elseif UE(i).ServingeNB == 3
        UE(i).Interference_dBm = UE(i).PowerReceive_dBm(1) + UE(i).PowerReceive_dBm(2);
        UE(i).PowerReceive_dBm = UE(i).PowerReceive_dBm(3);
    end
end

for i = 1:Num_UEs
    UE(i).SINR_dB = SINR_Calculation(UE(i).PowerReceive_dBm, UE(i).Interference_dBm, NoisePower_dBm);
end

%% CQI and Modulation Order Calculations

for i = 1:Num_UEs
	UE(i).CQI = CQI_Calculation(UE(i).SINR_dB);
    UE(i).Modulation = Modulation_Calculation(UE(i).CQI);
end

%% Assigning Services and Needed Resource Blocks to the UEs

for i = 1:Num_UEs
    UE(i).QCI = QCI_ExampleService();
    UE(i).ExampleService = ExampleService_Calculation(UE(i).QCI);
    [UE(i).Needed_RBs, UE(i).Needed_RBs_Type] = Needed_RBs_Calculation(UE(i).QCI);
end

%% Sorting UEs According to Distance

Distance1_3D = []; Allocated_UEs1 = 0; UE_Index1 = [];
Distance2_3D = []; Allocated_UEs2 = 0; UE_Index2 = [];
Distance3_3D = []; Allocated_UEs3 = 0; UE_Index3 = [];

for i = 1:Num_UEs
    if UE(i).ServingeNB == 1
        Allocated_UEs1 = Allocated_UEs1 + 1;
        UE_Index1(Allocated_UEs1, :) = i;
        Distance1_3D(Allocated_UEs1, :) = [UE(i).Distance(1), i];
    elseif UE(i).ServingeNB == 2
        Allocated_UEs2 = Allocated_UEs2 + 1;
        UE_Index2(Allocated_UEs2, :) = i;
        Distance2_3D(Allocated_UEs2, :) = [UE(i).Distance(2), i];
    elseif UE(i).ServingeNB == 3
        Allocated_UEs3 = Allocated_UEs3 + 1;
        UE_Index3(Allocated_UEs3, :) = i;
        Distance3_3D(Allocated_UEs3, :) = [UE(i).Distance(3), i];
    end
end

Sorted_Distance1_3D = sortrows(Distance1_3D, 'ascend');
For_Values1 = (Sorted_Distance1_3D(:, 2))';

Sorted_Distance2_3D = sortrows(Distance2_3D, 'ascend');
For_Values2 = (Sorted_Distance2_3D(:, 2))';

Sorted_Distance3_3D = sortrows(Distance3_3D, 'ascend');
For_Values3 = (Sorted_Distance3_3D(:, 2))';

%% Allocating Resource Blocks to the allocated UEs

Available_RBs1 = eNB(1).ResourceBlocks;  Limited_RBs1 = 0.3*Available_RBs1;
Available_RBs2 = eNB(2).ResourceBlocks;  Limited_RBs2 = 0.3*Available_RBs2;
Available_RBs3 = eNB(3).ResourceBlocks;  Limited_RBs3 = 0.3*Available_RBs3;

for i = For_Values1
    if UE(i).ServingeNB == 1
        if UE(i).Needed_RBs < Available_RBs1 && Available_RBs1 > Limited_RBs1
            UE(i).Allocated_RBs = UE(i).Needed_RBs;
            Available_RBs1 = Available_RBs1 - UE(i).Allocated_RBs;
        else
            UE(i).Allocated_RBs = 0;
        end
    end
end

for i = For_Values2
    if UE(i).ServingeNB == 2
        if UE(i).Needed_RBs < Available_RBs2 && Available_RBs2 > Limited_RBs2
            UE(i).Allocated_RBs = UE(i).Needed_RBs;
            Available_RBs2 = Available_RBs2 - UE(i).Allocated_RBs;
        else
            UE(i).Allocated_RBs = 0;
        end
    end
end

for i = For_Values3
    if UE(i).ServingeNB == 3
        if UE(i).Needed_RBs < Available_RBs3 && Available_RBs3 > Limited_RBs3
            UE(i).Allocated_RBs = UE(i).Needed_RBs;
            Available_RBs3 = Available_RBs3 - UE(i).Allocated_RBs;
        else
            UE(i).Allocated_RBs = 0;
        end
    end
end

%% Useful Calculations of the UEs

Location1 = []; RBs_BarChar1 = []; Allocated_UEs1 = 0; SINR1 = []; CQI1 = []; Mod1 = [];
Location2 = []; RBs_BarChar2 = []; Allocated_UEs2 = 0; SINR2 = []; CQI2 = []; Mod2 = [];
Location3 = []; RBs_BarChar3 = []; Allocated_UEs3 = 0; SINR3 = []; CQI3 = []; Mod3 = [];

for i = 1:Num_UEs
    if UE(i).ServingeNB == 1
        Allocated_UEs1 = Allocated_UEs1 + 1;
        Location1(Allocated_UEs1, :) = [Cord_x(i), Cord_y(i)];
        RBs_BarChar1(Allocated_UEs1, :) = [UE(i).Needed_RBs, UE(i).Allocated_RBs];
        SINR1(Allocated_UEs1, :) = UE(i).SINR_dB;
        CQI1(Allocated_UEs1, :) = UE(i).CQI;
        Mod1(Allocated_UEs1, :) = UE(i).Modulation;
    elseif UE(i).ServingeNB == 2
        Allocated_UEs2 = Allocated_UEs2 + 1;
        Location2(Allocated_UEs2, :) = [Cord_x(i), Cord_y(i)];
        RBs_BarChar2(Allocated_UEs2, :) = [UE(i).Needed_RBs, UE(i).Allocated_RBs];
        SINR2(Allocated_UEs2, :) = UE(i).SINR_dB;
        CQI2(Allocated_UEs2, :) = UE(i).CQI;
        Mod2(Allocated_UEs2, :) = UE(i).Modulation;
    elseif UE(i).ServingeNB == 3
        Allocated_UEs3 = Allocated_UEs3 + 1;
        Location3(Allocated_UEs3, :) = [Cord_x(i), Cord_y(i)];
        RBs_BarChar3(Allocated_UEs3, :) = [UE(i).Needed_RBs, UE(i).Allocated_RBs];
        SINR3(Allocated_UEs3, :) = UE(i).SINR_dB;
        CQI3(Allocated_UEs3, :) = UE(i).CQI;
        Mod3(Allocated_UEs3, :) = UE(i).Modulation;
    end
end

Distance1_Limit = []; Served_UEs1 = 0; Blocked_UEs1 = 0;
Distance2_Limit = []; Served_UEs2 = 0; Blocked_UEs2 = 0;
Distance3_Limit = []; Served_UEs3 = 0; Blocked_UEs3 = 0;

for i = For_Values1
    if UE(i).Allocated_RBs ~= 0
        Served_UEs1 = Served_UEs1 + 1;
        Distance1_Limit(Served_UEs1,:) = sqrt((eNB(1).x-UE(i).x)^2 + (eNB(1).y-UE(i).y)^2);

    else
        Blocked_UEs1 = Blocked_UEs1 + 1;
    end
end

for i = For_Values2
    if UE(i).Allocated_RBs ~= 0
        Served_UEs2 = Served_UEs2 + 1;
        Distance2_Limit(Served_UEs2,:) = sqrt((eNB(2).x-UE(i).x)^2 + (eNB(2).y-UE(i).y)^2);
    else
        Blocked_UEs2 = Blocked_UEs2 + 1;
    end
end

for i = For_Values3
    if UE(i).Allocated_RBs ~= 0
        Served_UEs3 = Served_UEs3 + 1;
        Distance3_Limit(Served_UEs3,:) = sqrt((eNB(3).x-UE(i).x)^2 + (eNB(3).y-UE(i).y)^2);
    else
        Blocked_UEs3 = Blocked_UEs3 + 1;
    end
end

theta = linspace(0, 2*pi, 100);

x1 = (eNB(1).x)./1e3 + 30*cos(theta);  y1 = (eNB(1).y)./1e3 + 25*sin(theta);
x2 = (eNB(2).x)./1e3 + 30*cos(theta);  y2 = (eNB(2).y)./1e3 + 25*sin(theta);
x3 = (eNB(3).x)./1e3 + 30*cos(theta);  y3 = (eNB(3).y)./1e3 + 25*sin(theta);

Max_Distance1_Limit = (max(Distance1_Limit))./1e3;
x4 = (eNB(1).x)./1e3 + Max_Distance1_Limit*cos(theta);  y4 = (eNB(1).y)./1e3 + Max_Distance1_Limit*sin(theta);

Max_Distance2_Limit = (max(Distance2_Limit))./1e3;
x5 = (eNB(2).x)./1e3 + Max_Distance2_Limit*cos(theta);  y5 = (eNB(2).y)./1e3 + Max_Distance2_Limit*sin(theta);

Max_Distance3_Limit = (max(Distance3_Limit))./1e3;
x6 = (eNB(3).x)./1e3 + Max_Distance3_Limit*cos(theta);  y6 = (eNB(3).y)./1e3 + Max_Distance3_Limit*sin(theta);

%% UEs vs Serving eNBs Plot

figure(1)

scatter(Cord_eNBx(1:3), Cord_eNBy(1:3), 300, 'filled', 'd', 'k'); hold on;
scatter(Location1(:,1), Location1(:,2), 100, 'filled', 'r', '<'); hold on;
scatter(Location2(:,1), Location2(:,2), 100, 'filled', 'b', '^'); hold on;
scatter(Location3(:,1), Location3(:,2), 100, 'filled', 'm', '>'); hold on;

plot(x1, y1); hold on; plot(x2, y2); hold on; plot(x3, y3); hold on;
plot(x5, y5, '--'); hold on; plot(x6, y6, '--'); hold on; plot(x4, y4, '--'); hold off;

grid on;
xlim([0 ((eNB(3).x)./1e3 + 5)]);  xlabel('Length (Km)');
ylim([0 ((eNB(3).y)./1e3 + 25)]);  ylabel('Width (Km)');
title('Plan View of eNBs and UEs');
legend('Serving eNBs', 'UEs served by eNB 1', 'UEs served by eNB 2', 'UEs served by eNB 3');

%% Resource Blocks Comparison Plot

figure(2)

subplot(3,1,1)
bar((UE_Index1)', RBs_BarChar1);
legend('Needed RBs','Allocated RBs'); grid on;
title('Resource Blocks Comparison for eNB 1');
xlabel('UEs Served by eNB 1'); ylabel('Resource Blocks');

subplot(3,1,2)
bar((UE_Index2)', RBs_BarChar2);
legend('Needed RBs','Allocated RBs'); grid on;
title('Resource Blocks Comparison for eNB 2');
xlabel('UEs Served by eNB 2'); ylabel('Resource Blocks');

subplot(3,1,3)
bar((UE_Index3)', RBs_BarChar3);
legend('Needed RBs','Allocated RBs'); grid on;
title('Resource Blocks Comparison for eNB 3');
xlabel('UEs Served by eNB 3'); ylabel('Resource Blocks');

%% SINR Plot

figure(3)

subplot(3,1,1)
bar((UE_Index1)', SINR1, 'FaceColor', [0.494 0.184 0.556]);
title('SINR of UEs Served by eNB 1'); grid on;
xlabel('UEs Served by eNB 1'); ylabel('SINR (dB)');

subplot(3,1,2)
bar((UE_Index2)', SINR2, 'FaceColor', [0.494 0.184 0.556]);
title('SINR of UEs Served by eNB 2'); grid on;
xlabel('UEs Served by eNB 2'); ylabel('SINR (dB)');

subplot(3,1,3)
bar((UE_Index3)', SINR3, 'FaceColor', [0.494 0.184 0.556]);
title('SINR of UEs Served by eNB 3'); grid on;
xlabel('UEs Served by eNB 3'); ylabel('SINR (dB)');

%% CQI Plot

figure(4)

subplot(3,1,1)
bar((UE_Index1)', CQI1, 'FaceColor', [0.850 0.325 0.098]);
title('CQI of UEs Served by eNB 1'); grid on;
xlabel('UEs Served by eNB 1'); ylabel('CQI');

subplot(3,1,2)
bar((UE_Index2)', CQI2, 'FaceColor', [0.850 0.325 0.098]);
title('CQI of UEs Served by eNB 2'); grid on;
xlabel('UEs Served by eNB 2'); ylabel('CQI');

subplot(3,1,3)
bar((UE_Index3)', CQI3, 'FaceColor', [0.850 0.325 0.098]);
title('CQI of UEs Served by eNB 3'); grid on;
xlabel('UEs Served by eNB 3'); ylabel('CQI');

%% Modulation Order Plot

figure(5)

subplot(3,1,1)
bar((UE_Index1)', Mod1, 'FaceColor', [0.466 0.674 0.188]);
title('Modulation Order of UEs Served by eNB 1'); grid on;
xlabel('UEs Served by eNB 1'); ylabel('Modulation Oreder');

subplot(3,1,2)
bar((UE_Index2)', Mod2, 'FaceColor', [0.466 0.674 0.188]);
title('Modulation Order of UEs Served by eNB 2'); grid on;
xlabel('UEs Served by eNB 2'); ylabel('Modulation Oreder');

subplot(3,1,3)
bar((UE_Index3)', Mod3, 'FaceColor', [0.466 0.674 0.188]);
title('Modulation Order of UEs Served by eNB 3'); grid on;
xlabel('UEs Served by eNB 3'); ylabel('Modulation Oreder');
%}
%}