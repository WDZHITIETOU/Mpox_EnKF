classdef MyStochasticEnKF < handle
% MyStochasticEnKF The Ensemble Kalman Filter whose update step will be
    % carried out stochastically.
    % Detailed instructions are shown here.

    properties
    % all properties in MyStochasticEnKF class
        epsilon = 1e-1 .* ones(1,4); % Lower bound of infection rate
        N0 (1,1) double {mustBePositive} = 1e8 % Total number of people simulated
        stepSize (1,1) double {mustBePositive} = 1 % Time step predicted by the model
        Ts (1,1) double {mustBePositive} = 1 % Time step for data correction
        omega (1,1) double {mustBePositive} = 1/7.19 % Reciprocal of the incubation period
        gamma (1,1) double {mustBePositive} = 1/7 % Reciprocal of the infectious period
        nState (1,1) int8 {mustBePositive} = 32 % Number of states
        nObservation (1,1) int8 {mustBePositive} = 4 % The number of observations

        mu0 (1,32) double % The initial mean of the state vector
        Sigma0 (32,32) double % The initial covariance of the state vector
        sampleSize (1,1) int16 % Number of samples

        MeasurementNoiseCovariance0 (4,4) double % Initial
        ProcessNoiseCovariance0 (32,32) double % Initial 
        MeasurementNoiseCovariance (4,4) double % Time-varying
        ProcessNoiseCovariance (32,32) double % Time-varying


        x (32,:) double % The state matrix at each time step
        x0 (32,:) double % The initial state matrix at each time step
        originalData
        observations  
        filteredEstimations

        t % Days of observation
        d datetime % Date of observation
        tt % Days of status update
        dd datetime % Date of status update

        % Recurrence Curves
        time
        SS
        EE
        II
        RR
        IncidenceRate
    end
    
    methods
        function md = MyStochasticEnKF()
            %MYSTOCHASTICENKF constructor
            % Detailed instructions are shown here

            %import observed data
            caseDataTable = readtable('../Raw data/assimilatedData.xlsx');
            %summary(caseDataTable);
            caseDataTable.Properties.VariableNames(2) = "date";
            caseDataTable.age_group = categorical(caseDataTable.age_group);
            caseDataTable = removevars(caseDataTable, "time");

            uniqueDates = sort(unique(caseDataTable.date));
            allDates = min(uniqueDates):max(uniqueDates);
            intersection = setxor(allDates,uniqueDates);

            caseData1 = unstack(caseDataTable,'I','age_group','VariableNamingRule','preserve');
            addData = cell(numel(intersection),width(caseData1));
            intersection = num2cell(intersection); 
            addData(:,1) = intersection;
            addData(:,2:end) = {0};
            %[addData{:,2:end}] = deal(0);
            caseData1 = [caseData1;addData];
            caseData1 = sortrows(caseData1,1);
            caseData1 = fillmissing(caseData1,'constant',0,'DataVariables',[2 3 4 5]);
            
            uniqueDates1 = sort(unique(caseData1.date));
            d0 = min(uniqueDates1);
            md.t = days(uniqueDates1 - d0);
            md.d = d0 + days(md.t);

            md.originalData = caseData1{:,[2 3 4 5]};
            md.observations = caseData1{:,[2 3 4 5]};
            md.observations = smoothdata(md.observations, 1, 'gaussian', 5);
            
            % Setting initial values for initializing sampleSize samples
            S0 = [md.N0*0.3021, md.N0*0.3932, md.N0*0.2085, md.N0*0.0963];
            E0 = [1, 20, 10, 1]; % mustBeNonZeros in updateMeasurementFcn
            I0 = [1, 20, 10, 1];
            Beta = 1e-1 * ones(4);
            Beta([1,end], :) = Beta([1,end], :) * 1e-1;
            md.mu0 = [S0, ... % S
                      E0,... % E
                      I0,... % I
                      zeros(1,4),...% R
                      Beta(:)']; % Beta
            md.Sigma0 = diag([1e6*ones(1,4), repmat([1 9 4 1],1,2), ones(1,4),... % S E I R
                             1e-3*ones(1,4), 1e-2*ones(1,8), 1e-3*ones(1,4)]); % Beta

            s = 1e1*ones(1,24); % The variance of each component of the derivatives
            s([6 7 10 11 18 19 22 23]) = s([6 7 10 11 18 19 22 23]) * 1e1;
            Q = zeros(md.nState-16); % [S(1-4), E(1-4), I(1-4), R(1-4)]
            Q(1,5) = -sum(s(1:4));
            Q(2,6) = -sum(s(5:8));
            Q(3,7) = -sum(s(9:12));
            Q(4,8) = -sum(s(13:16));
            Q(5,9) = -s(17);
            Q(6,10) = -s(18);
            Q(7,11) = -s(19);
            Q(8,12) = -s(20);
            Q(9,13) = -s(21);
            Q(10,14) = -s(22);
            Q(11,15) = -s(23);
            Q(12,16) = -s(24);

            Q = Q + diag([sum(s(1:4)), sum(s(5:8)), sum(s(9:12)), sum(s(13:16)),...
                sum(s(1:4))+s(17), sum(s(5:8))+s(18), sum(s(9:12))+s(19), sum(s(13:16))+s(20),...
                s(17)+s(21), s(18)+s(22), s(19)+s(23), s(20)+s(24),...
                s(21), s(22), s(23), s(24)]) / 2;
            Q = Q + Q.';
            Q(Q < 1e-8) = 0;
            Q = blkdiag(Q, 1e-6*eye(16));

            md.ProcessNoiseCovariance0 = Q;
            md.ProcessNoiseCovariance = md.ProcessNoiseCovariance0;
            md.MeasurementNoiseCovariance0 = diag([4; 1e2; 25; 1])*1e-1;
            md.MeasurementNoiseCovariance = md.MeasurementNoiseCovariance0;
        end

        function filtering(md)
            n = md.Ts/md.stepSize;
            md.x0 = initializeSamples(md);
            md.x = md.x0;

            for k = 1:size(md.observations,1)
                DataCorrection(md, md.observations(k,:).');
                
                for j = 1:n
                    md.filteredEstimations((k-1)*n+j,:) = mean(md.x, 2);
                    ModelPrediction(md);
                    UpdateProcessNoiseCovariance(md);
                end
                UpdateMeasurementNoiseCovariance(md);
            end

            md.tt = md.stepSize * (0 : n*size(md.observations,1)-1);
            md.dd = md.d(1) + days(md.tt - md.tt(1));
           
        end

        function UpdateMeasurementNoiseCovariance(md)
            incidenceRate = mean(MeasurementFcn(md),2);
            E0 = md.mu0(5:8);
            amplifier = log( incidenceRate(:) ./ E0(:) ); 
            md.MeasurementNoiseCovariance = diag(amplifier) * md.MeasurementNoiseCovariance0 * diag(amplifier)';
        end

        function UpdateProcessNoiseCovariance(md)       
            dxdt0 = abs(mean(computeDerivatives(md, md.x0), 2));
            dxdt = abs(mean(computeDerivatives(md, md.x), 2));
            idx = 1:16;
            amplifier = log(dxdt(idx) ./ dxdt0(idx));
            md.ProcessNoiseCovariance(idx,idx) = diag(amplifier) * md.ProcessNoiseCovariance0(idx,idx) * diag(amplifier)';
            % if any(isnan(md.ProcessNoiseCovariance))
            %     1;
            % end
        end

        function ModelPrediction(md)
            w = mvnrnd(zeros(md.nState,1), md.ProcessNoiseCovariance, md.sampleSize)';
            md.x = StateTransitionFcn(md) + w;
            %md.x(17:end,:) = md.x(17:end,:) * 0.99;
            md.x = max(0, md.x);
        end

        function DataCorrection(md, y)
            v = mvnrnd(zeros(md.nObservation,1), md.MeasurementNoiseCovariance, md.sampleSize)';
            S = cov(md.x');
            C = S * 1.2; % variance inflation
            %C = C; % localization
            H = MeasurementJacobian(md);
            gain = C*H'/(H*C*H' + md.MeasurementNoiseCovariance);
            md.x = md.x +  gain*(y + v - MeasurementFcn(md));
            md.x = max(0, md.x);
        end
         
        function x0 = initializeSamples(md)
            x0 = mvnrnd(md.mu0, md.Sigma0, md.sampleSize)';
            x0 = max(0,x0);       
        end

        function y = MeasurementFcn(md)
            E = md.x(5:8,:);
            y = md.omega.*E;
        end

        function H = MeasurementJacobian(md)
            H = zeros(4, size(md.x,1));
            H(1,5) = md.omega;
            H(2,6) = md.omega;
            H(3,7) = md.omega;
            H(4,8) = md.omega;
        end

        function xNext = StateTransitionFcn(md)
            dxdt = computeDerivatives(md, md.x);
            xNext = md.x + md.stepSize * dxdt; 
            xNext = max(0, xNext);
        end

        function dxdt = computeDerivatives(md,x)
            % compute derivatives: Calculate the first derivative of each state over time
            % Detailed instructions are shown here
            % sampleSize samples are restored in an n*r matrix x
            S1 = x(1,:);
            S2 = x(2,:);
            S3 = x(3,:);
            S4 = x(4,:);
            E1 = x(5,:);
            E2 = x(6,:);
            E3 = x(7,:);
            E4 = x(8,:);
            I1 = x(9,:);
            I2 = x(10,:);
            I3 = x(11,:);
            I4 = x(12,:);
            R1 = x(13,:);
            R2 = x(14,:);
            R3 = x(15,:);
            R4 = x(16,:);
            N1 = sum(x([1 5 9 13],:));
            N2 = sum(x([2 6 10 14],:));
            N3 = sum(x([3 7 11 15],:));
            N4 = sum(x([4 8 12 16],:));
            beta11 = x(17,:);
            beta12 = x(18,:);
            beta13 = x(19,:);
            beta14 = x(20,:);
            beta21 = x(21,:);
            beta22 = x(22,:);
            beta23 = x(23,:);
            beta24 = x(24,:);
            beta31 = x(25,:);
            beta32 = x(26,:);
            beta33 = x(27,:);
            beta34 = x(28,:);
            beta41 = x(29,:);
            beta42 = x(30,:);
            beta43 = x(31,:);
            beta44 = x(32,:);

            dxdt = zeros(size(x));
            dxdt(1,:) = -beta11.*S1.*I1./N1 - beta21.*S1.*I2./N1 - beta31.*S1.*I3./N1 - beta41.*S1.*I4./N1 - md.epsilon(1);
            dxdt(2,:) = -beta12.*S2.*I1./N2 - beta22.*S2.*I2./N2 - beta32.*S2.*I3./N2 - beta42.*S2.*I4./N2 - md.epsilon(2);
            dxdt(3,:) = -beta13.*S3.*I1./N3 - beta23.*S3.*I2./N3 - beta33.*S3.*I3./N3 - beta43.*S3.*I4./N3 - md.epsilon(3);
            dxdt(4,:) = -beta14.*S4.*I1./N4 - beta24.*S4.*I2./N4 - beta34.*S4.*I3./N4 - beta44.*S4.*I4./N4 - md.epsilon(4);
            dxdt(5,:) = beta11.*S1.*I1./N1 + beta21.*S1.*I2./N1 + beta31.*S1.*I3./N1 + beta41.*S1.*I4./N1 + md.epsilon(1) - md.omega.*E1;
            dxdt(6,:) = beta12.*S2.*I1./N2 + beta22.*S2.*I2./N2 + beta32.*S2.*I3./N2 + beta42.*S2.*I4./N2 + md.epsilon(2) - md.omega.*E2;
            dxdt(7,:) = beta13.*S3.*I1./N3 + beta23.*S3.*I2./N3 + beta33.*S3.*I3./N3 + beta43.*S3.*I4./N3 + md.epsilon(3) - md.omega.*E3;
            dxdt(8,:) = beta14.*S4.*I1./N4 + beta24.*S4.*I2./N4 + beta34.*S4.*I3./N4 + beta44.*S4.*I4./N4 + md.epsilon(4) - md.omega.*E4;
            dxdt(9,:) = md.omega.*E1 - md.gamma.*I1;
            dxdt(10,:) = md.omega.*E2 - md.gamma.*I2;
            dxdt(11,:) = md.omega.*E3 - md.gamma.*I3;
            dxdt(12,:) = md.omega.*E4 - md.gamma.*I4;
            dxdt(13,:) = md.gamma.*I1;
            dxdt(14,:) = md.gamma.*I2;
            dxdt(15,:) = md.gamma.*I3;
            dxdt(16,:) = md.gamma.*I4;
        end

        function simulateWithFilteredBeta(md, tGrid, BetaMatrix, initialStates, tSpan)
            % INPUT:
            %   md - model object
            %   initialStates - [S(:); E(:); I(:); R(:)]
            %   tSpan - 1*2 array, [tStart, tEnd]
            BetaFunction = griddedInterpolant(tGrid, BetaMatrix, 'linear', 'previous');  
            S.type = '()';
            %S.subs = 1:16;
            S.subs = {1:16,':'}; 
            odefun = @(t, y) subsref(md.computeDerivatives([y; BetaFunction(t)']), S);
            [md.time, y] = ode45(odefun, tSpan, initialStates);
            md.SS = y(:,1:4);
            md.EE = y(:,5:8);
            md.II = y(:,9:12);
            md.RR = y(:,13:16);
            md.IncidenceRate = md.omega * y(:,5:8);
        end

        function R = computeReff(md, beta11, beta12, beta13, beta14, beta21, beta22, beta23, beta24, beta31, beta32, beta33, beta34, beta41, beta42, beta43, beta44)
            % computeReff: Calculate the effective reproduction number of each age group
            % Detailed instructions are shown here
            R = zeros(16+1, size(beta11,1));
            R(1,:) = beta11 ./ md.gamma;
            R(2,:) = beta12 ./ md.gamma; 
            R(3,:) = beta13 ./ md.gamma;
            R(4,:) = beta14 ./ md.gamma;
            R(5,:) = beta21 ./ md.gamma;
            R(6,:) = beta22 ./ md.gamma;
            R(7,:) = beta23 ./ md.gamma;
            R(8,:) = beta24 ./ md.gamma;
            R(9,:) = beta31 ./ md.gamma;
            R(10,:) = beta32 ./ md.gamma;
            R(11,:) = beta33 ./ md.gamma;
            R(12,:) = beta34 ./ md.gamma;
            R(13,:) = beta41 ./ md.gamma;
            R(14,:) = beta42 ./ md.gamma;
            R(15,:) = beta43 ./ md.gamma;
            R(16,:) = beta44 ./ md.gamma;
            proportionGroup1 = 0.3021;
            proportionGroup2 = 0.3932;
            proportionGroup3 = 0.2085;
            proportionGroup4 = 0.0963;
            R(17,:) = (R(1,:)+R(2,:)+R(3,:)+R(4,:)).*proportionGroup1 + (R(5,:)+R(6,:)+R(7,:)+R(8,:)).*proportionGroup2 + (R(9,:)+R(10,:)+R(11,:)+R(12,:)).*proportionGroup3 + (R(13,:)+R(14,:)+R(15,:)+R(16,:)).*proportionGroup4;
        end
        
        function visualize(md)

            % observation

            fig1 = figure("Name", "observation");
            fig1.WindowState = 'maximized';
            fig1.Color = [1,1,1];
            an = tiledlayout(4,2);
            an.YLabel.String = 'Daily Incidence (individuals per day)';
            an.YLabel.FontSize = 24;
            an.XLabel.String = 'Year/Month';
            an.XLabel.FontSize = 24;
            an.YLabel.FontWeight = 'bold';
            an.XLabel.FontWeight = 'bold';
            an.XLabel.FontName = 'Times New Roman';
            an.YLabel.FontName = 'Times New Roman';

            ax = nexttile(2);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,1), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(4);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,2), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(6);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,3), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(8);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,4), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(1);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.originalData(:,1), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(3);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';            
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.originalData(:,2), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(5);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.originalData(:,3), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile(7);
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.originalData(:,4), "Color", "#0072BD", "LineWidth", 2);
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off
            
            % beta(t)
            %Different colors represent susceptible individuals of different ages
            %Different linestyles represent infectious individuals of different ages

            fig2 = figure("Color", [1,1,1], "Name", "beta(t)");
            fig2.WindowState = "maximized";
            T = tiledlayout(2,2,Padding="compact");
            T.YLabel.String = '\beta_{ij}';
            T.YLabel.FontSize = 24;
            T.XLabel.String = 'Year/Month';
            T.XLabel.FontSize = 24;
            T.YLabel.FontWeight = 'bold';
            T.XLabel.FontWeight = 'bold';
            T.XLabel.FontName = "Times New Roman";
            T.YLabel.FontName = "Times New Roman";

            for i = 1:4
                ax = nexttile;
                ax.FontName = "Times New Roman";
                ax.FontWeight = "bold";
                ax.FontSize = 16;
                ax.LineWidth = 1;
                ax.Box = "on";
                hold on
                plot(md.d, md.filteredEstimations(:, [16+i, 20+i, 24+i, 28+i]), LineWidth=2);
                lgd = legend("\beta_{1" + i +"}", "\beta_{2" + i + "}", "\beta_{3" + i + "}", "\beta_{4" + i + "}");
                lgd.Box = "off";
                lgd.FontSize = 14;
                ax.XAxis.TickLabelFormat = "u/M";
                hold off
            end                                                                               

            % prediction and observation

            fig3 = figure('Name', 'prediction and observation');
            fig3.WindowState = 'maximized';
            fig3.Color = [1,1,1];
            an = tiledlayout(4,1);
            an.YLabel.String = 'Daily Incidence (individuals per day)';
            an.YLabel.FontSize = 24;
            an.XLabel.String = 'Year/Month';
            an.XLabel.FontSize = 24;
            an.YLabel.FontWeight = 'bold';
            an.XLabel.FontWeight = 'bold';
            an.YLabel.FontName = 'Times New Roman';
            an.XLabel.FontName = 'Times New Roman';

            ax = nexttile;
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.filteredEstimations(:,5)*md.omega, "LineWidth", 3);
            sca = scatter(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,1), 80, "x", LineWidth=1.5);
            sca.MarkerEdgeAlpha = 0.7;
            lgd = legend('observation', 'prediction');
            lgd.Box = "off";
            lgd.FontSize = 20;
            lgd.FontWeight = "bold";
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile;
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.filteredEstimations(:,6)*md.omega, "LineWidth", 3);
            sca = scatter(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,2), 80, "x", LineWidth=1.5);
            sca.MarkerEdgeAlpha = 0.7;
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile;
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.filteredEstimations(:,7)*md.omega, "LineWidth", 3);
            sca = scatter(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,3), 80, "x", LineWidth=1.5);
            sca.MarkerEdgeAlpha = 0.7;
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off

            ax = nexttile;
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.TickLength = [0.005,0.005];
            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.FontName = 'Times New Roman';
            hold on
            plot(datetime(2022,4,12):datetime(2023,2,15), md.filteredEstimations(:,8)*md.omega, "LineWidth", 3);
            sca = scatter(datetime(2022,4,12):datetime(2023,2,15), md.observations(:,4), 80, "x", LineWidth=1.5);
            sca.MarkerEdgeAlpha = 0.7;
            ax.XAxis.TickLabelFormat = 'u/M';
            hold off
            
            % R(t)
            %Different colors represent susceptible individuals of different ages
            %Different linestyles represent infectious individuals of different ages
          
            Rt = computeReff(md,md.filteredEstimations(:,17),md.filteredEstimations(:,18),md.filteredEstimations(:,19),md.filteredEstimations(:,20),...
                md.filteredEstimations(:,21),md.filteredEstimations(:,22),md.filteredEstimations(:,23),md.filteredEstimations(:,24),...
                md.filteredEstimations(:,25),md.filteredEstimations(:,26),md.filteredEstimations(:,27),md.filteredEstimations(:,28),...
                md.filteredEstimations(:,29),md.filteredEstimations(:,30),md.filteredEstimations(:,31),md.filteredEstimations(:,32));
            
            xvalues = {'0-17','18-44','45-64','65+'};
            yvalues = {'0-17','18-44','45-64','65+'};
            fig4 = figure;
            fig4.WindowState = 'maximized';
            fig4.Color = [1,1,1];
            T = tiledlayout(3,3,Padding="compact");
            T.YLabel.String = 'To age group j';
            T.YLabel.FontSize = 24;
            T.YLabel.FontName = "Times New Roman";
            T.XLabel.String = 'From age group i';
            T.XLabel.FontSize = 24;
            T.XLabel.FontName = "Times New Roman";
            T.YLabel.FontWeight = 'bold';
            T.XLabel.FontWeight = 'bold';
            
            riqi = ["2022/5/15", "2022/6/15", "2022/7/15", "2022/8/15", "2022/9/15",...
                    "2022/10/15", "2022/11/15", "2022/12/15", "2023/1/15"];
            index = [34, 65, 95, 126, 157, 187, 218, 248, 279];
            for i = 1:9
                nexttile;
                Rmatrix = reshape(Rt(1:16,index(i)),[4 4]);
                h = heatmap(xvalues,yvalues,Rmatrix);
                h.FontName = "Times New Roman";
                h.FontSize = 16;
                h.Title = riqi(i);
            end
                                                                                        
            fig5 = figure(Color=[1 1 1], Name="R(t)");
            fig5.WindowState = 'maximized';
            T = tiledlayout(2,2);
            T.YLabel.String = 'R_{ij}';
            T.YLabel.FontSize = 24;
            T.XLabel.String = 'Year/Month';
            T.XLabel.FontSize = 24;
            T.YLabel.FontWeight = 'bold';
            T.XLabel.FontWeight = 'bold';
            T.XLabel.FontName = "Times New Roman";
            T.YLabel.FontName = "Times New Roman";

            for i = 1:4
                ax = nexttile;
                ax.FontName = "Times New Roman";
                ax.FontWeight = "bold";
                ax.FontSize = 16;
                ax.LineWidth = 1;
                ax.Box = "on";
                hold on
                plot(md.d, Rt([i, 4+i, 8+i, 12+i],:), LineWidth=2);
                lgd = legend("R_{1" + i +"}", "R_{2" + i + "}", "R_{3" + i + "}", "R_{4" + i + "}");
                lgd.Box = "off";
                lgd.FontSize = 14;
                ax.XAxis.TickLabelFormat = "u/M";
                hold off
            end
        end
    end
end

