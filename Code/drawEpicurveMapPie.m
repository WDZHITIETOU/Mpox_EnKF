clc;
clear;
close all;

%读取世界国家行政区划矢量文件
countries = readgeotable("../Raw data/map/map/World_countries.shp");

countries.FCNAME(220) = "刚果";
countries.FCNAME(226) = "刚果";

%读取世界各国猴痘病例确诊数
opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Country", "WHORegion", "TotalConfirmedCases", "TotalProbableCases", "TotalDeaths"];
opts.VariableTypes = ["string", "categorical", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "skip";
opts = setvaropts(opts, "Country", "WhitespaceRule", "trim");
opts = setvaropts(opts, ["Country", "WHORegion"], "EmptyFieldRule", "auto");
mpoxCasesbyCountry = readtable("../Raw data/mpox cases and deaths by country.csv", opts);

countries = outerjoin(countries, mpoxCasesbyCountry, "LeftKeys", "FCNAME", "RightKeys", "Country", "MergeKeys", true, "RightVariables", "TotalConfirmedCases", "Type", "left");
countries.TotalConfirmedCases = fillmissing(countries.TotalConfirmedCases, "constant", 0);
    
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
tiledlayout(5,5,Padding="compact");
newcolors = ["#F0EEEF", "#C6CCDC", "#9DACCB",...
             "#7789B7", "#CBD7C3", "#ACBF9F",...
             "#89AA7B", "#D6D6D6", "#EB6969"];
colororder(F, newcolors);

nexttile(1, [3 4]);
mx = newmap;
geoplot(countries, "ColorVariable", "TotalConfirmedCases")
mx.LineWidth = 1;
mx.FontName = "Times New Roman";
mx.FontSize = 12;
mx.FontWeight = "bold";
colormap(mx, sky)
cb = colorbar;
cb.Label.String = "Total Confirmed Cases";
cb.FontName = "Times New Roman";
cb.FontSize = 12;
cb.FontWeight = 'bold';

nexttile(5);
[~, idx] = max(mpoxCasesbyRegionandAgeGroup{1, 2:end});
explode = false(1,9);
explode(idx) = true;
pie1 = piechart(mpoxCasesbyRegionandAgeGroup{1, 2:end}, ExplodedWedges=explode);
pie1.LabelStyle = "none";
pie1.FaceAlpha = 1;
pie1.LegendVisible = "on";
pie1.LegendTitle = "Age Group";
names = {'0-9','10-17','18-29','30-39','40-49','50-59','60-69','70-79','80+'};
pie1.Names = names;
pie1.FontSize = 14;
pie1.FontName = "Times New Roman";

nexttile(10)
[~, idx] = max(mpoxCasesbyRegionandAgeGroup{2, 2:end});
explode = false(1,9);
explode(idx) = true;
pie2 = piechart(mpoxCasesbyRegionandAgeGroup{2, 2:end}, ExplodedWedges=explode);
pie2.LabelStyle = "none";
pie2.FaceAlpha = 1;

nexttile(15)
[~, idx] = max(mpoxCasesbyRegionandAgeGroup{3, 2:end});
explode = false(1,9);
explode(idx) = true;
pie3 = piechart(mpoxCasesbyRegionandAgeGroup{3, 2:end}, ExplodedWedges=explode);
pie3.LabelStyle = "none";
pie3.FaceAlpha = 1;

bar1 = nexttile(16, [2 4]);
bar(dailyIncidencebyAgeGroup.yearWeek, dailyIncidencebyAgeGroup{:, [2 3 4 5]}, 1,'stacked');
xticks(["2022-1" "2022-10" "2022-20" "2022-30" "2022-40" "2023-1"]);
xlabel("Year-Week")
ylabel("Incidence (individuals per day)")
legend('0-17', '18-44', '45-64', '65+')
newcolors2 = [140 191 135;
              62 96 141;
              203 148 117;
              144 146 145] ./ 255;
colororder(bar1, newcolors2);
bar1.FontName = "Times New Roman";
bar1.FontSize = 12;
bar1.FontWeight = "bold";
bar1.Box = "off";
bar1.LineWidth = 1;

nexttile(20)
[~, idx] = max(mpoxCasesbyRegionandAgeGroup{4, 2:end});
explode = false(1,9);
explode(idx) = true;
pie4 = piechart(mpoxCasesbyRegionandAgeGroup{4, 2:end}, ExplodedWedges=explode);
pie4.LabelStyle = "none";
pie4.FaceAlpha = 1;

nexttile(25)
[~, idx] = max(mpoxCasesbyRegionandAgeGroup{5, 2:end});
explode = false(1,9);
explode(idx) = true;
pie5 = piechart(mpoxCasesbyRegionandAgeGroup{5, 2:end}, ExplodedWedges=explode);
pie5.LabelStyle = "none";
pie5.FaceAlpha = 1;

%%
% 创建 textbox
annotation(F,'textbox',...
    [0.11172362685266 0.932692307692307 0.0251551874455099 0.0705128205128203],...
    'String',{'a'},...
    'Margin',4,...
    'LineWidth',1,...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');

% 创建 textbox
annotation(F,'textbox',...
    [0.112595466434176 0.443910256410256 0.0251551874455099 0.076923076923077],...
    'String','b',...
    'Margin',4,...
    'LineWidth',1,...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');

% 创建 textbox
annotation(F,'textbox',...
    [0.791758500435921 0.894230769230769 0.0242833478639927 0.076923076923077],...
    'String','c',...
    'Margin',4,...
    'LineWidth',1,...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');