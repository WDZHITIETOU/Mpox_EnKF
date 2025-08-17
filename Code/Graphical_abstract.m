clc
clear
close all

%读取全球猴痘病例每日新增数
opts = spreadsheetImportOptions("NumVariables", 4);
opts.Sheet = "要拟合的数据CSV";
opts.DataRange = "A2:D995";
opts.VariableNames = ["ageGroup", "datetime", "days", "incidence"];
opts.VariableTypes = ["categorical", "datetime", "double", "double"];
opts = setvaropts(opts, "ageGroup", "EmptyFieldRule", "auto");
dailyIncidencebyAgeGroup = readtable("../Raw data/daily incidence by age group.xlsx", opts, "UseExcel", false);

dailyIncidencebyAgeGroup.year = year(dailyIncidencebyAgeGroup.datetime);
dailyIncidencebyAgeGroup.week = week(dailyIncidencebyAgeGroup.datetime);
yearWeek = cell(height(dailyIncidencebyAgeGroup), 1);
for i = 1:height(dailyIncidencebyAgeGroup)
    yearWeek{i} = append(num2str(dailyIncidencebyAgeGroup.year(i)), '-', num2str(dailyIncidencebyAgeGroup.week(i)));
end
valueset = ["2022-" + string(1:53), "2023-" + string(1:7)];
dailyIncidencebyAgeGroup.yearWeek = categorical(yearWeek, valueset);
dailyIncidencebyAgeGroup = pivot(dailyIncidencebyAgeGroup, Rows="yearWeek", Columns="ageGroup", DataVariable="incidence", Method="sum", IncludeEmptyGroups=true);

load ../'Intermediate data'/estimate_beta.mat
meanR = zeros(17,310);
lowerLimitR = zeros(17,310);
upperLimitR = zeros(17,310);
for i = 1:17
    for j = 1:310
        temp = reshape(R(i,j,:), [], 1);
        pd = fitdist(temp, "Normal");
        meanR(i,j) = mean(pd);
        lowerLimitR(i,j) = icdf(pd, 0.025);
        upperLimitR(i,j) = icdf(pd, 0.975);
    end
end

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

figure = figure('Color',[1 1 1]);
figure.WindowState = 'maximized';
tiledlayout(1,3,Padding="compact");

%ax = nexttile(1, [3 3]);
ax = nexttile;
ax.FontName = "Times New Roman";
ax.FontSize = 14;
ax.FontWeight = "bold";
ax.Box = "on";
ax.LineWidth = 1;
hold on
bar(dailyIncidencebyAgeGroup.yearWeek, dailyIncidencebyAgeGroup{:, [2 3 4 5]}, 1,'stacked');
xticks(["2022-1" "2022-10" "2022-20" "2022-30" "2022-40" "2023-1"]);
xlabel("Year-Week")
ylabel("Incidence (individuals per day)")
legend('0-17', '18-44', '45-64', '65+')
newcolors2 = [140 191 135;
              62 96 141;
              203 148 117;
              144 146 145] ./ 255;
colororder(ax, newcolors2);
hold off

%ax = nexttile(5, [3 3]);
ax = nexttile;
ax.FontName = "Times New Roman";
ax.FontSize = 14;
ax.FontWeight = "bold";
ax.Box = "on";
ax.LineWidth = 1;
hold on
plot(datetime(2022,4,12):datetime(2023,2,15), meanR(6,:), "LineWidth", 2);
ax.XAxis.TickLabelFormat = 'M/u';
xlabel("Month/Year")
ylabel("Effective reproduction number")
xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
fill(xdata,[lowerLimitR(6,:), flip(upperLimitR(6,:))],[0 0.4470 0.7410],...
            FaceAlpha = 0.3, LineStyle = "none");
R1 = yline(1,"--","R = 1","LineWidth", 1, "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment", "center");
R1.FontName = "Times New Roman";
R1.FontSize = 14;
R1.FontWeight = "bold";
R2 = xline(datetime(2022,5,1), "--", "burn-in", "LineWidth", 1, "LabelVerticalAlignment", "middle", "LabelHorizontalAlignment", "left", "LabelOrientation", "aligned");
R2.FontName = "Times New Roman";
R2.FontSize = 14;
R2.FontWeight = "bold";
hold off

%ax = nexttile(9, [3 3]);
ax = nexttile;
ax.FontName = "Times New Roman";
ax.FontSize = 14;
ax.FontWeight = "bold";
ax.Box = "on";
ax.LineWidth = 1;
xticks(ax, 1:6)
xticklabels(ax, ["90%" "80%" "60%" "40%" "20%" "no intervention"])
hold on
errorbar(1:5, vaccination_mean(1,:,1), vaccination_neg(1,:,1), vaccination_pos(1,:,1), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
errorbar(1:5, vaccination_mean(1,:,2), vaccination_neg(1,:,2), vaccination_pos(1,:,2), 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
errorbar(6, noIntervention_mean, noIntervention_neg, noIntervention_pos, 'o', 'LineWidth', 2, 'CapSize', 20, 'MarkerSize', 10)
xlabel("Vaccination coverage")
ylabel("Cumulative case count")
hold off