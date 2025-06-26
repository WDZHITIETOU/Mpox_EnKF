clear
clc
close all

infectiousPeriodSample = [0 4 7 10 28];
cdInfectiousPeriod = [0, 0.25, 0.5, 0.75, 1];
pdInfectious = makedist('PiecewiseLinear','x',infectiousPeriodSample,'Fx',cdInfectiousPeriod);

incubationPeriodSample1 = [8.5*ones(1,18), 7.6*ones(1,40), 7.0*ones(1,23), 6.0*ones(1,77), 7.0*ones(1,144), 9.1*ones(1,30), 6.0*ones(1,112), 11.0*ones(1,16),...
                          13.0*ones(1,14), 7.0*ones(1,527), 8.0*ones(1,78), 5.6*ones(1,35), 7.8*ones(1,54), 7.0*ones(1,29), 8.2*ones(1,209)]';
pdIncubation = fitdist(incubationPeriodSample1, 'Kernel', 'Support', 'positive');

pdN0 = makedist('Uniform', 'lower', 1e8, 'upper', 2e8);

load ../'Intermediate data'/globalStream.mat
RandStream.setGlobalStream(globalStream)
myStream = RandStream.getGlobalStream;
myStream.State = myState;
lhs = lhsdesign(100,3);

filteredX = zeros(310,32,100);
filteredRT = zeros(17,310,100);
for i = 1:100
    RandStream.setGlobalStream(globalStream)
    myStream = RandStream.getGlobalStream;
    myStream.State = myState;

    md = MyStochasticEnKF;
    md.sampleSize = 1000;
    md.epsilon = [1, 1e2, 1e2, 1] * 1e-1;
    md.N0 = floor(icdf(pdN0, lhs(i,1)));
    md.omega = 1 ./ icdf(pdIncubation, lhs(i,2));
    md.gamma = 1 ./ icdf(pdInfectious, lhs(i,3));
    filtering(md);

    Rt = computeReff(md,md.filteredEstimations(:,17),md.filteredEstimations(:,18),md.filteredEstimations(:,19),md.filteredEstimations(:,20),...
        md.filteredEstimations(:,21),md.filteredEstimations(:,22),md.filteredEstimations(:,23),md.filteredEstimations(:,24),...
        md.filteredEstimations(:,25),md.filteredEstimations(:,26),md.filteredEstimations(:,27),md.filteredEstimations(:,28),...
        md.filteredEstimations(:,29),md.filteredEstimations(:,30),md.filteredEstimations(:,31),md.filteredEstimations(:,32));
    filteredX(:,:,i) = md.filteredEstimations;
    filteredRT(:,:,i) = Rt;

end

save ../'Intermediate data'/sensitivity_analysis.mat filteredX filteredRT

%%
clc
clear
close all

load ../'Intermediate data'/sensitivity_analysis.mat

median_filteredX = median(filteredX, 3);
%max_filteredX = max(filteredX, [], 3);
%min_filteredX = min(filteredX, [], 3);
max_filteredX = quantile(filteredX, 0.75, 3);
min_filteredX = quantile(filteredX, 0.25, 3);

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
         4 8 12 16] + 16;
name = ["\beta_{11}", "\beta_{21}", "\beta_{31}", "\beta_{41}";...
        "\beta_{12}", "\beta_{22}", "\beta_{32}", "\beta_{42}";...
        "\beta_{13}", "\beta_{23}", "\beta_{33}", "\beta_{43}";...
        "\beta_{14}", "\beta_{24}", "\beta_{34}", "\beta_{44}"];

for i = 1:4
    ax = nexttile;
    ax.FontName = "Times New Roman";
    ax.FontWeight = "bold";
    ax.FontSize = 14;
    ax.Box = "on";
    ax.LineWidth = 1;
    ax.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
    hold on
    plot(datetime(2022,4,12):datetime(2023,2,15), median_filteredX(:, index(i,:)), LineWidth = 1.5);
    %ylim()
    ax.XAxis.TickLabelFormat = 'u/M';
    lgd = legend(name(i,:));
    lgd.AutoUpdate = "off";
    lgd.Box = "off";
    lgd.FontSize = 14;
    lgd.FontWeight = "bold";
    lgd.FontName = "Times New Roman";
    xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
    fill(xdata,[min_filteredX(:,index(i,1)); flip(max_filteredX(:,index(i,1)))],[0 0.4470 0.7410],...
         xdata,[min_filteredX(:,index(i,2)); flip(max_filteredX(:,index(i,2)))],[0.8500 0.3250 0.0980],...
         xdata,[min_filteredX(:,index(i,3)); flip(max_filteredX(:,index(i,3)))],[0.9290 0.6940 0.1250],...
         xdata,[min_filteredX(:,index(i,4)); flip(max_filteredX(:,index(i,4)))],[0.4940 0.1840 0.5560],...
         FaceAlpha = 0.3, LineStyle = "none");
end

median_filteredRT = median(filteredRT, 3);
%max_filteredRT = max(filteredRT, [], 3);
%min_filteredRT = min(filteredRT, [], 3);
max_filteredRT = quantile(filteredRT, 0.75, 3);
min_filteredRT = quantile(filteredRT, 0.25, 3);

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
        ax.FontSize = 14;
        ax.Box = "on";
        ax.LineWidth = 1;
        ax.ColorOrder = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
        hold on
        plot(datetime(2022,4,12):datetime(2023,2,15), median_filteredRT([index(2.*i-1) index(2.*i)],:), "LineWidth", 1.5);
        ylim([0 2])
        ax.XAxis.TickLabelFormat = 'u/M';
        lgd = legend(name(2.*i-1), name(2.*i));
        lgd.AutoUpdate = "off";
        lgd.Box = "off";
        lgd.FontSize = 14;
        lgd.FontWeight = "bold";
        lgd.FontName = "Times New Roman";
        xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
        fill(xdata,[min_filteredRT(index(i.*2-1),:), flip(max_filteredRT(index(i.*2-1),:))],[0 0.4470 0.7410],...
             xdata,[min_filteredRT(index(i.*2),:), flip(max_filteredRT(index(i.*2),:))],[0.8500 0.3250 0.0980],...
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
    else
        ax = nexttile;
        ax.FontName = "Times New Roman";
        ax.FontWeight = "bold";
        ax.FontSize = 14;
        ax.Box = "on";
        ax.LineWidth = 1;
        hold on
        plot(datetime(2022,4,12):datetime(2023,2,15), median_filteredRT(17,:), "LineWidth", 1.5);
        ylim([0 2])
        ax.XAxis.TickLabelFormat = 'u/M';
        lgd = legend("R_{eff}");
        lgd.AutoUpdate = "off";
        lgd.Box = "off";
        lgd.FontSize = 14;
        lgd.FontWeight = "bold";
        lgd.FontName = "Times New Roman";
        xdata = [datetime(2022,4,12):datetime(2023,2,15), datetime(2023,2,15):-1:datetime(2022,4,12)];
        fill(xdata,[min_filteredRT(17,:), flip(max_filteredRT(17,:))],[0 0.4470 0.7410],...
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
    end
end
