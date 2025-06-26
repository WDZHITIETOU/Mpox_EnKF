classdef intervention < handle
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Beta
        omega = 1/7.19
        timePoint = 1
        gamma
        gamma1 = 1/7
        gamma2 = 1/7
        daysCount = 310-50
    end
    
    methods
        function obj = intervention(beta)
            %UNTITLED 构造此类的实例
            %   此处显示详细说明
            obj.Beta = beta(51:end,:);
        end
        
        function dydt = odefun(obj,t,y)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            %   y是列向量
            if t < obj.timePoint
                obj.gamma = obj.gamma1;
            elseif t >= obj.timePoint && t <= 92
                obj.gamma = obj.gamma2;
            else
                obj.gamma = obj.gamma1;
            end

            S = y(1:4);
            E = y(5:8);
            I = y(9:12);
            R = y(13:16);
            N = [sum(y([1 5 9 13])); sum(y([2 6 10 14])); sum(y([3 7 11 15])); sum(y([4 8 12 16]))]; 
            BetaMatrix = reshape(interp1((1:obj.daysCount)', obj.Beta, t), 4, 4);
            infectionRate = (BetaMatrix .* (S./N)) * I;
            incidenceRate = obj.omega .* E;
            removeRate = obj.gamma .* I;
            dSdt = -infectionRate;
            dEdt = infectionRate - incidenceRate;
            dIdt = incidenceRate - removeRate;
            dRdt = removeRate;
            dydt = [dSdt;dEdt;dIdt;dRdt];
        end
        
        function [cumulativeCases, peakTimeofNewCases] = computeIndices(obj, X)
            % description
            y0 = X(51,1:16);
            tSpan = [1, obj.daysCount];
            [t,y] = ode45(@(t,y) odefun(obj,t,y), tSpan, y0);
            tt = (1:obj.daysCount)';
            y = interp1(t,y,tt);

            cumulativeCases = sum(obj.omega.*y(:,5:8),'all') + sum(y0(9:16));                                        
            [~, peakTimeofNewCases] = max(sum(obj.omega.*y(:,5:8),2));

        end
    end
end