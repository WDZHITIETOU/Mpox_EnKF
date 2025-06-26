clear;
clc;
close all;

% globalStream = RandStream.getGlobalStream;
% myState = globalStream.State;
% save ../'Intermediate data'/globalStream.mat

load ../'Intermediate data'/globalStream.mat
RandStream.setGlobalStream(globalStream)
myStream = RandStream.getGlobalStream;
myStream.State = myState;

md = MyStochasticEnKF;
md.omega = 1./7.19;
md.gamma = 1./7;
md.N0 = 1e8;
md.sampleSize = 1000;
md.epsilon = [1, 1e2, 1e2, 1] .* 1e-1;
filtering(md);
visualize(md);
