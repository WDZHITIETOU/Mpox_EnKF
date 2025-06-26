clear
clc
close all

load ../'Intermediate data'/globalStream.mat
RandStream.setGlobalStream(globalStream)
myStream = RandStream.getGlobalStream;
myStream.State = myState;

X = zeros(310,16,100);
beta = zeros(310,16,100);
R = zeros(17,310,100);
for i = 1:100
    md = MyStochasticEnKF;
    md.omega = 1./7.19;
    md.gamma = 1./7;
    md.N0 = 1e8;
    md.sampleSize = 1000;
    md.epsilon = [1, 1e2, 1e2, 1] .* 1e-1;
    filtering(md);
    Rt = computeReff(md,md.filteredEstimations(:,17),md.filteredEstimations(:,18),md.filteredEstimations(:,19),md.filteredEstimations(:,20),...
        md.filteredEstimations(:,21),md.filteredEstimations(:,22),md.filteredEstimations(:,23),md.filteredEstimations(:,24),...
        md.filteredEstimations(:,25),md.filteredEstimations(:,26),md.filteredEstimations(:,27),md.filteredEstimations(:,28),...
        md.filteredEstimations(:,29),md.filteredEstimations(:,30),md.filteredEstimations(:,31),md.filteredEstimations(:,32));
    X(:,:,i) = md.filteredEstimations(:,1:16);
    beta(:,:,i) = md.filteredEstimations(:,17:end);
    R(:,:,i) = Rt; 
end

save ../'Intermediate data'/estimate_beta.mat X beta R

%%
clc
clear
close all

load ../'Intermediate data'/estimate_beta.mat

meanBeta = zeros(310,16);
lowerLimitBeta = zeros(310,16);
upperLimitBeta = zeros(310,16);
for i = 1:310
    for j = 1:16
        temp = reshape(beta(i,j,:), [], 1);
        pd = fitdist(temp, "Normal");
        meanBeta(i,j) = mean(pd);
        lowerLimitBeta(i,j) = icdf(pd, 0.025);
        upperLimitBeta(i,j) = icdf(pd, 0.975);
    end
end

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

figure1 = figure('Color',[1 1 1]);
figure1.WindowState = 'maximized';
T = tiledlayout(2,2,Padding="compact");
T.YLabel.String = '\beta_{ij}';
T.YLabel.FontSize = 30;
T.YLabel.FontName = "Times New Roman";
T.XLabel.String = 'Year/Month';
T.XLabel.FontSize = 30;
T.XLabel.FontName = "Times New Roman";
T.YLabel.FontWeight = 'bold';
T.XLabel.FontWeight = 'bold';

index = [1 5 9 13;...
         2 6 10 14;...
         3 7 11 15;...
         4 8 12 16];
name = ["\beta_{11}", "\beta_{21}", "\beta_{31}", "\beta_{41}";...
        "\beta_{12}", "\beta_{22}", "\beta_{32}", "\beta_{42}";...
        "\beta_{13}", "\beta_{23}", "\beta_{33}", "\beta_{43}";...
        "\beta_{14}", "\beta_{24}", "\beta_{34}", "\beta_{44}"];

for i = 1:4
    ax = nexttile;
    ax.FontName = "Times New Roman";
    ax.FontWeight = "bold";
    ax.FontSize = 18;
    ax.Box = "on";
    ax.LineWidth = 1;
    ax.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
    hold on
    plot(datetime(2022,4,12):datetime(2023,2,15), meanBeta(:, index(i,:)), LineWidth = 2);
    %ylim()
    ax.XAxis.TickLabelFormat = 'u/M';
    lgd = legend(name(i,:));
    lgd.AutoUpdate = "off";
    lgd.Box = "off";
    lgd.FontSize = 18;
    lgd.FontWeight = "bold";
    lgd.FontName = "Times New Roman";
    xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
    fill(xdata,[lowerLimitBeta(:,index(i,1)); flip(upperLimitBeta(:,index(i,1)))],[0 0.4470 0.7410],...
         xdata,[lowerLimitBeta(:,index(i,2)); flip(upperLimitBeta(:,index(i,2)))],[0.8500 0.3250 0.0980],...
         xdata,[lowerLimitBeta(:,index(i,3)); flip(upperLimitBeta(:,index(i,3)))],[0.9290 0.6940 0.1250],...
         xdata,[lowerLimitBeta(:,index(i,4)); flip(upperLimitBeta(:,index(i,4)))],[0.4940 0.1840 0.5560],...
         FaceAlpha = 0.3, LineStyle = "none");
end

figure2 = figure('Color',[1 1 1]);
figure2.WindowState = 'maximized';
T = tiledlayout(3,3,Padding="compact");
T.YLabel.String = 'R_{ij}';
T.YLabel.FontSize = 30;
T.YLabel.FontName = "Times New Roman";
T.XLabel.String = 'Year/Month';
T.XLabel.FontSize = 30;
T.XLabel.FontName = "Times New Roman";
T.YLabel.FontWeight = 'bold';
T.XLabel.FontWeight = 'bold';

index = [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16];
name = ["R_{11}", "R_{21}", "R_{31}", "R_{41}",...
        "R_{12}", "R_{22}", "R_{32}", "R_{42}",...
        "R_{13}", "R_{23}", "R_{33}", "R_{43}",...
        "R_{14}", "R_{24}", "R_{34}", "R_{44}"];
for i = 1:9
    if i ~= 9
        ax = nexttile;
        ax.FontName = "Times New Roman";
        ax.FontWeight = "bold";
        ax.FontSize = 18;
        ax.Box = "on";
        ax.LineWidth = 1;
        ax.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
        hold on
        plot(datetime(2022,4,12):datetime(2023,2,15), meanR([index(2.*i-1) index(2.*i)],:), "LineWidth", 2);
        ylim([0 2])
        ax.XAxis.TickLabelFormat = 'u/M';
        lgd = legend(name(2.*i-1), name(2.*i));
        lgd.AutoUpdate = "off";
        lgd.Box = "off";
        lgd.FontSize = 18;
        lgd.FontWeight = "bold";
        lgd.FontName = "Times New Roman";
        xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
        fill(xdata,[lowerLimitR(index(i.*2-1),:), flip(upperLimitR(index(i.*2-1),:))],[0 0.4470 0.7410],...
             xdata,[lowerLimitR(index(i.*2),:), flip(upperLimitR(index(i.*2),:))],[0.8500 0.3250 0.0980],...
             FaceAlpha = 0.3, LineStyle = "none");
        R1 = yline(1,"--","R = 1","LineWidth", 1, "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment", "center");
        R1.FontName = "Times New Roman";
        R1.FontSize = 18;
        R1.FontWeight = "bold";
        R2 = xline(datetime(2022,5,1), "--", "burn-in", "LineWidth", 1, "LabelVerticalAlignment", "middle", "LabelHorizontalAlignment", "left", "LabelOrientation", "aligned");
        R2.FontName = "Times New Roman";
        R2.FontSize = 18;
        R2.FontWeight = "bold";
        hold off
    else
        ax = nexttile;
        ax.FontName = "Times New Roman";
        ax.FontWeight = "bold";
        ax.FontSize = 18;
        ax.Box = "on";
        ax.LineWidth = 1;
        hold on
        plot(datetime(2022,4,12):datetime(2023,2,15), meanR(17,:), "LineWidth", 2);
        ylim([0 2])
        ax.XAxis.TickLabelFormat = 'u/M';
        lgd = legend("R_{eff}");
        lgd.AutoUpdate = "off";
        lgd.Box = "off";
        lgd.FontSize = 18;
        lgd.FontWeight = "bold";
        lgd.FontName = "Times New Roman";
        xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
        fill(xdata,[lowerLimitR(17,:), flip(upperLimitR(17,:))],[0 0.4470 0.7410],...
             FaceAlpha = 0.3, LineStyle = "none");
        R1 = yline(1,"--","R = 1","LineWidth", 1, "LabelVerticalAlignment", "bottom", "LabelHorizontalAlignment", "center");
        R1.FontName = "Times New Roman";
        R1.FontSize = 18;
        R1.FontWeight = "bold";
        R2 = xline(datetime(2022,5,1), "--", "burn-in", "LineWidth", 1, "LabelVerticalAlignment", "middle", "LabelHorizontalAlignment", "left", "LabelOrientation", "aligned");
        R2.FontName = "Times New Roman";
        R2.FontSize = 18;
        R2.FontWeight = "bold";
        hold off
    end
end

xvalues = {'0-17','18-44','45-64','65+'};
yvalues = {'0-17','18-44','45-64','65+'};
fig4 = figure;
fig4.WindowState = 'maximized';
fig4.Color = [1,1,1];
T = tiledlayout(3,3,Padding="compact");
T.YLabel.String = 'To age group j';
T.YLabel.FontSize = 30;
T.YLabel.FontName = "Times New Roman";
T.XLabel.String = 'From age group i';
T.XLabel.FontSize = 30;
T.XLabel.FontName = "Times New Roman";
T.YLabel.FontWeight = 'bold';
T.XLabel.FontWeight = 'bold';

riqi = ["2022/5/15", "2022/6/15", "2022/7/15", "2022/8/15", "2022/9/15",...
    "2022/10/15", "2022/11/15", "2022/12/15", "2023/1/15"];
index = [34, 65, 95, 126, 157, 187, 218, 248, 279];
for i = 1:9
    nexttile;
    Rmatrix = reshape(meanR(1:16,index(i)),[4 4]);
    h = heatmap(xvalues,yvalues,Rmatrix);
    h.FontName = "Times New Roman";
    h.FontSize = 18;
    h.Title = riqi(i);
end