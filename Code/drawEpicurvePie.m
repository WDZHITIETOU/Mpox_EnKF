clc;
clear;
close all;

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

%读取世界各大区分年龄分性别的猴痘病例数
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["sex", "ageGroup", "region", "cases"];
opts.VariableTypes = ["categorical", "categorical", "categorical", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "skip";
opts = setvaropts(opts, ["sex", "ageGroup", "region"], "EmptyFieldRule", "auto");
mpoxCasesbyRegionandAgeGroup = readtable("../Raw data/mpox cases age-sex pyramid by WHO region.csv", opts);

mpoxCasesbyRegionandAgeGroup = pivot(mpoxCasesbyRegionandAgeGroup, Rows="region", Columns="ageGroup", DataVariable="cases", Method="sum");

%绘图
F = figure('Color',[1 1 1]);
F.WindowState = "maximized";
tiledlayout(5,6,Padding="compact");
newcolors = ["#F0EEEF", "#C6CCDC", "#9DACCB",...
             "#7789B7", "#CBD7C3", "#ACBF9F",...
             "#89AA7B", "#D6D6D6", "#EB6969"];
colororder(F, newcolors);

ax = nexttile(1, [5 5]);
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
ax.FontName = "Times New Roman";
ax.FontSize = 14;
ax.FontWeight = "bold";
ax.Box = "off";
ax.LineWidth = 1;

for i = 1:5
    ax1 = nexttile(6*i);
    [~, idx] = max(mpoxCasesbyRegionandAgeGroup{i, 2:end});
    explode = false(1,9);
    explode(idx) = true;
    pie = piechart(mpoxCasesbyRegionandAgeGroup{i, 2:end}, ExplodedWedges=explode);
    pie.LabelStyle = "none";
    pie.FaceAlpha = 1;
    pie.FontSize = 14;
    pie.FontName = "Times New Roman";
    if i == 1
        pie.LegendVisible = "on";
        pie.LegendTitle = "Age Group";
        names = {'0-9','10-17','18-29','30-39','40-49','50-59','60-69','70-79','80+'};        
        pie.Names = names;
    end
end
