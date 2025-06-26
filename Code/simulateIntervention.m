clc
clear
close all

load ../'Intermediate data'/estimate_beta.mat

noIntervention_cumulativeCases = zeros(1,100);
noIntervention_peakTime = zeros(1,100);

vaccination_whole_cumulativeCases = zeros(6,5,100);
vaccination_highRisk_cumulativeCases = zeros(6,5,100);
vaccination_whole_peakTime = zeros(6,5,100);
vaccination_highRisk_peakTime = zeros(6,5,100);

interruptTransmission_whole_cumulativeCases = zeros(6,5,100);
interruptTransmission_highRisk_cumulativeCases = zeros(6,5,100);
interruptTransmission_whole_peakTime = zeros(6,5,100);
interruptTransmission_highRisk_peakTime = zeros(6,5,100);

caseManagement_whole_cumulativeCases = zeros(6,6,100);
caseManagement_highRisk_cumulativeCases = zeros(6,6,100);
caseManagement_whole_peakTime = zeros(6,6,100);
caseManagement_highRisk_peakTime = zeros(6,6,100);

timeIndex = {1:92, 15:92, 31:92, 45:92, 62:92, 76:92};
coverageIndex = [0.9, 0.8, 0.6, 0.4, 0.2];
reductionIndex = [0.9, 0.8, 0.6, 0.4, 0.2];
timePointIndex = [1, 15, 31, 45, 62, 76];
gammaIndex = [1, 1/2, 1/3, 1/4, 1/5, 1/6];

% noIntervention
for i = 1:100
    obj = intervention(beta(:,:,i));
    [noIntervention_cumulativeCases(i), noIntervention_peakTime(i)] = computeIndices(obj, X(:,:,i));
end

% vaccination
for i = 1:100
    for j = 1:6
        for k = 1:5
            obj = intervention(beta(:,:,i));
            obj.Beta(timeIndex{j}, :) = obj.Beta(timeIndex{j}, :) .* (1 - 0.8.*coverageIndex(k));
            [vaccination_whole_cumulativeCases(j,k,i), vaccination_whole_peakTime(j,k,i)] = computeIndices(obj, X(:,:,i));
        end
    end
end
for i = 1:100
    for j = 1:6
        for k = 1:5
            obj = intervention(beta(:,:,i));
            obj.Beta(timeIndex{j}, [2 6 10 14]) = obj.Beta(timeIndex{j}, [2 6 10 14]) .* (1 - 0.8.*coverageIndex(k));
            [vaccination_highRisk_cumulativeCases(j,k,i), vaccination_highRisk_peakTime(j,k,i)] = computeIndices(obj, X(:,:,i));
        end
    end
end

% interruptTransmission
for i = 1:100
    for j = 1:6
        for k = 1:5
            obj = intervention(beta(:,:,i));
            obj.Beta(timeIndex{j}, :) = obj.Beta(timeIndex{j}, :) .* (1 - reductionIndex(k));
            [interruptTransmission_whole_cumulativeCases(j,k,i), interruptTransmission_whole_peakTime(j,k,i)] = computeIndices(obj, X(:,:,i));
        end
    end
end
for i = 1:100
    for j = 1:6
        for k = 1:5
            obj = intervention(beta(:,:,i));
            obj.Beta(timeIndex{j}, [2 5 6 7 8 10 14]) = obj.Beta(timeIndex{j}, [2 5 6 7 8 10 14]) .* (1 - reductionIndex(k));
            [interruptTransmission_highRisk_cumulativeCases(j,k,i), interruptTransmission_highRisk_peakTime(j,k,i)] = computeIndices(obj, X(:,:,i));
        end
    end
end

% caseManagement
for i = 1:100
    for j = 1:6
        for k = 1:6
            obj = intervention(beta(:,:,i));
            obj.timePoint = timePointIndex(j);
            obj.gamma2 = gammaIndex(k);
            [caseManagement_whole_cumulativeCases(j,k,i), caseManagement_whole_peakTime(j,k,i)] = computeIndices(obj, X(:,:,i));
        end
    end
end
for i = 1:100
    for j = 1:6
        for k = 1:6
            obj = intervention(beta(:,:,i));
            obj.timePoint = timePointIndex(j);
            obj.gamma2 = [1/7; gammaIndex(k); 1/7; 1/7];
            [caseManagement_highRisk_cumulativeCases(j,k,i), caseManagement_highRisk_peakTime(j,k,i)] = computeIndices(obj, X(:,:,i));
        end
    end
end

save ../'Intermediate data'/simulate_intervention.mat noIntervention_cumulativeCases noIntervention_peakTime ...
    vaccination_whole_cumulativeCases vaccination_whole_peakTime vaccination_highRisk_cumulativeCases ...
    vaccination_highRisk_peakTime interruptTransmission_whole_cumulativeCases interruptTransmission_whole_peakTime ...
    interruptTransmission_highRisk_cumulativeCases interruptTransmission_highRisk_peakTime caseManagement_whole_cumulativeCases ...
    caseManagement_whole_peakTime caseManagement_highRisk_cumulativeCases caseManagement_highRisk_peakTime

%%
clc
clear
close all

load ../'Intermediate data'/simulate_intervention.mat

% no intervention
pd = fitdist(noIntervention_cumulativeCases', "Normal");
noIntervention_mean = mean(pd);
if isnan(icdf(pd, 0.025))
    noIntervention_neg = 0;
else
    noIntervention_neg = mean(pd) - icdf(pd, 0.025);
end
if isnan(icdf(pd, 0.975))
    noIntervention_pos = 0;
else
    noIntervention_pos = icdf(pd, 0.975) - mean(pd);
end

% vaccination cumulativeCases
vaccination_mean = zeros(6,5,2);
vaccination_neg = zeros(6,5,2);
vaccination_pos = zeros(6,5,2);
for i = 1:6
    for j = 1:5
        temp = reshape(vaccination_whole_cumulativeCases(i,j,:),[],1);
        pd = fitdist(temp, "Normal");
        vaccination_mean(i,j,1) = mean(pd);
        if isnan(icdf(pd, 0.025))
            vaccination_neg(i,j,1) = 0;
        else
            vaccination_neg(i,j,1) = mean(pd) - icdf(pd, 0.025);
        end
        if isnan(icdf(pd, 0.975))
            vaccination_pos(i,j,1) = 0;
        else
            vaccination_pos(i,j,1) = icdf(pd, 0.975) - mean(pd);
        end
    end
end
for i = 1:6
    for j = 1:5
        temp = reshape(vaccination_highRisk_cumulativeCases(i,j,:),[],1);
        pd = fitdist(temp, "Normal");
        vaccination_mean(i,j,2) = mean(pd);
        if isnan(icdf(pd, 0.025))
            vaccination_neg(i,j,2) = 0;
        else
            vaccination_neg(i,j,2) = mean(pd) - icdf(pd, 0.025);
        end
        if isnan(icdf(pd, 0.975))
            vaccination_pos(i,j,2) = 0;
        else
            vaccination_pos(i,j,2) = icdf(pd, 0.975) - mean(pd);
        end
    end
end

% interruptTransmission cumulativeCases
interruptTransmission_mean = zeros(6,5,2);
interruptTransmission_neg = zeros(6,5,2);
interruptTransmission_pos = zeros(6,5,2);
for i = 1:6
    for j = 1:5
        temp = reshape(interruptTransmission_whole_cumulativeCases(i,j,:),[],1);
        pd = fitdist(temp, "Normal");
        interruptTransmission_mean(i,j,1) = mean(pd);
        if isnan(icdf(pd, 0.025))
            interruptTransmission_neg(i,j,1) = 0;
        else
            interruptTransmission_neg(i,j,1) = mean(pd) - icdf(pd, 0.025);
        end
        if isnan(icdf(pd, 0.975))
            interruptTransmission_pos(i,j,1) = 0;
        else
            interruptTransmission_pos(i,j,1) = icdf(pd, 0.975) - mean(pd);
        end
    end
end
for i = 1:6
    for j = 1:5
        temp = reshape(interruptTransmission_highRisk_cumulativeCases(i,j,:),[],1);
        pd = fitdist(temp, "Normal");
        interruptTransmission_mean(i,j,2) = mean(pd);
        if isnan(icdf(pd, 0.025))
            interruptTransmission_neg(i,j,2) = 0;
        else
            interruptTransmission_neg(i,j,2) = mean(pd) - icdf(pd, 0.025);
        end
        if isnan(icdf(pd, 0.975))
            interruptTransmission_pos(i,j,2) = 0;
        else
            interruptTransmission_pos(i,j,2) = icdf(pd, 0.975) - mean(pd);
        end
    end
end

% caseManagement
caseManagement_mean = zeros(6,6,2);
caseManagement_neg = zeros(6,6,2);
caseManagement_pos = zeros(6,6,2);
for i = 1:6
    for j = 1:6
        temp = reshape(caseManagement_whole_cumulativeCases(i,j,:),[],1);
        pd = fitdist(temp, "Normal");
        caseManagement_mean(i,j,1) = mean(pd);
        if isnan(icdf(pd, 0.025))
            caseManagement_neg(i,j,1) = 0;
        else
            caseManagement_neg(i,j,1) = mean(pd) - icdf(pd, 0.025);
        end
        if isnan(icdf(pd, 0.975))
            caseManagement_pos(i,j,1) = 0;
        else
            caseManagement_pos(i,j,1) = icdf(pd, 0.975) - mean(pd);
        end
    end
end
for i = 1:6
    for j = 1:6
        temp = reshape(caseManagement_highRisk_cumulativeCases(i,j,:),[],1);
        pd = fitdist(temp, "Normal");
        caseManagement_mean(i,j,2) = mean(pd);
        if isnan(icdf(pd, 0.025))
            caseManagement_neg(i,j,2) = 0;
        else
            caseManagement_neg(i,j,2) = mean(pd) - icdf(pd, 0.025);
        end
        if isnan(icdf(pd, 0.975))
            caseManagement_pos(i,j,2) = 0;
        else
            caseManagement_pos(i,j,2) = icdf(pd, 0.975) - mean(pd);
        end
    end
end

% plot
figure1 = figure('Color',[1 1 1]);
figure1.WindowState = 'maximized';
T = tiledlayout(2,3,Padding="compact");
T.YLabel.String = 'Cumulative case count';
T.YLabel.FontSize = 30;
T.YLabel.FontName = "Times New Roman";
T.XLabel.String = "Vaccination coverage";
T.XLabel.FontSize = 30;
T.XLabel.FontName = "Times New Roman";
T.YLabel.FontWeight = 'bold';
T.XLabel.FontWeight = 'bold';
for i = 1:6
    ax = nexttile;
    ax.FontName = "Times New Roman";
    ax.FontWeight = "bold";
    ax.FontSize = 18;
    ax.Box = "on";
    ax.LineWidth = 1;
    xlim(ax, [0.5 6.5])
    ylim(ax, [0 7e4])
    xticks(ax, 1:6)
    xticklabels(ax, ["90%" "80%" "60%" "40%" "20%" "no intervention"])
    hold on
    errorbar(1:5, vaccination_mean(i,:,1), vaccination_neg(i,:,1), vaccination_pos(i,:,1), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    errorbar(1:5, vaccination_mean(i,:,2), vaccination_neg(i,:,2), vaccination_pos(i,:,2), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    errorbar(6, noIntervention_mean, noIntervention_neg, noIntervention_pos, 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    hold off
end
lgd = legend(["Intervention targeting whole population", "Intervention targeting high-risk groups", "No intervention"]);
lgd.FontSize = 18;
lgd.FontWeight = "bold";
lgd.FontName = "Times New Roman";
lgd.Orientation = "horizontal";
lgd.Layout.Tile = "south";


figure2 = figure('Color',[1 1 1]);
figure2.WindowState = 'maximized';
T = tiledlayout(2,3,Padding="compact");
T.YLabel.String = 'Cumulative case count';
T.YLabel.FontSize = 30;
T.YLabel.FontName = "Times New Roman";
T.XLabel.String = "Percentage reduction in transmission";
T.XLabel.FontSize = 30;
T.XLabel.FontName = "Times New Roman";
T.YLabel.FontWeight = 'bold';
T.XLabel.FontWeight = 'bold';
for i = 1:6
    ax = nexttile;
    ax.FontName = "Times New Roman";
    ax.FontWeight = "bold";
    ax.FontSize = 18;
    ax.Box = "on";
    ax.LineWidth = 1;
    xlim(ax, [0.5 6.5])
    ylim(ax, [0 7e4])
    xticks(ax, 1:6)
    xticklabels(ax, ["90%" "80%" "60%" "40%" "20%" "no intervention"])
    hold on
    errorbar(1:5, interruptTransmission_mean(i,:,1), interruptTransmission_neg(i,:,1), interruptTransmission_pos(i,:,1), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    errorbar(1:5, interruptTransmission_mean(i,:,2), interruptTransmission_neg(i,:,2), interruptTransmission_pos(i,:,2), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    errorbar(6, noIntervention_mean, noIntervention_neg, noIntervention_pos, 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    hold off
end
lgd = legend(["Intervention targeting whole population", "Intervention targeting high-risk groups", "No intervention"]);
lgd.FontSize = 18;
lgd.FontWeight = "bold";
lgd.FontName = "Times New Roman";
lgd.Orientation = "horizontal";
lgd.Layout.Tile = "south";

figure3 = figure('Color',[1 1 1]);
figure3.WindowState = 'maximized';
T = tiledlayout(2,3,Padding="compact");
T.YLabel.String = 'Cumulative case count';
T.YLabel.FontSize = 30;
T.YLabel.FontName = "Times New Roman";
T.XLabel.String = "The infectious period is shorten to";
T.XLabel.FontSize = 30;
T.XLabel.FontName = "Times New Roman";
T.YLabel.FontWeight = 'bold';
T.XLabel.FontWeight = 'bold';
for i = 1:6
    ax = nexttile;
    ax.FontName = "Times New Roman";
    ax.FontWeight = "bold";
    ax.FontSize = 18;
    ax.Box = "on";
    ax.LineWidth = 1;
    xlim(ax, [0.5 7.5])
    ylim(ax, [0 7e4])
    xticks(ax, 1:7)
    xticklabels(ax, ["1 day" "2 days" "3 days" "4 days" "5 days" "6 days" "no intervention"])
    hold on
    errorbar(1:6, caseManagement_mean(i,:,1), caseManagement_neg(i,:,1), caseManagement_pos(i,:,1), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    errorbar(1:6, caseManagement_mean(i,:,2), caseManagement_neg(i,:,2), caseManagement_pos(i,:,2), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    errorbar(7, noIntervention_mean, noIntervention_neg, noIntervention_pos, 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
    hold off
end
lgd = legend(["Intervention targeting whole population", "Intervention targeting high-risk groups", "No intervention"]);
lgd.FontSize = 18;
lgd.FontWeight = "bold";
lgd.FontName = "Times New Roman";
lgd.Orientation = "horizontal";
lgd.Layout.Tile = "south";
