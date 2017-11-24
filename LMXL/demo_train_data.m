
clear all
clc

load('videodata100.mat');

INPUT.Xa = XMAT(:,4:end);
EstimOpt.NamesA = {'price';'commer';'fast';'more TV';'more both';'share use'; 'share all';'no service'};

INPUT.Y = XMAT(:,3);

EstimOpt.NP = 1e2;
EstimOpt.NAlt = 5;
EstimOpt.NCT = 11;
EstimOpt.NRep = 2000; %This is what Train uses
EstimOpt.Draws = 6;
EstimOpt.NObs = EstimOpt.NP*EstimOpt.NCT;

% EstimOpt.Dist = zeros(size(INPUT.Xa,2),1);
% EstimOpt.Dist(end) = 1;

% EstimOpt.WTP_space = 0;
% INPUT.Xa = [INPUT.Xa(:,2:end),-INPUT.Xa(:,1)/10];
% EstimOpt.NamesA = {'commer';'fast';'more TV';'more both';'share use'; 'share all';'no service';'-price/10'};

EstimOpt.Order = 3; % this is what Train uses by default
EstimOpt.Dist = 2;

% This is what Train uses - mean +/- 2 s.d. from the MXL_d(?) model
% EstimOpt.Bounds = [Results.MXL_d.bhat(1:8) - 2*Results.MXL_d.bhat(9:16).^2.^0.5,Results.MXL_d.bhat(1:8) + 2*Results.MXL_d.bhat(9:16).^2.^0.5];
EstimOpt.Bounds = [-0.804619871971509,0.459131971110239;-0.205483280809719,0.114142064395781;-1.63710844563277,1.03743894056728;-2.29577994925115,1.48395847930476;-1.04513037070863,1.81347030453542;-1.74989971176418,1.66666419885748;-0.688765208764180,1.91391608771255;-11.7910451803257,0.683978272739412];

[INPUT,Results,EstimOpt,OptimOpt] = DataCleanDCE(INPUT,EstimOpt);

EstimOpt.NumGrad = 0;
Results.MNL = MNL(INPUT,Results,EstimOpt,OptimOpt);
% EstimOpt.FullCov = 0;
% Results.MXL_d = MXL(INPUT,Results,EstimOpt,OptimOpt);
% EstimOpt.FullCov = 1;
% Results.MXL = MXL(INPUT,Results,EstimOpt,OptimOpt);

EstimOpt.NumGrad = 1;
EstimOpt.FullCov =0;
Results.LMXL_d = LMXL(INPUT,Results,EstimOpt,OptimOpt);
EstimOpt.FullCov =1;
Results.LMXL = LMXL(INPUT,Results,EstimOpt,OptimOpt);
