function [histwgt,midbin] = histwgt(vv, ww, startpt, endpt, nbins) 

%Inputs: 
%vv: column vector of values 
%ww: column vestor of weights 
%startpt - lowest value to use for histogram 
%endpt - highest value to use for histogram.
%nbins - number of bins 
%Histogram has nbins evenly spaced bins between startpt and endpt. 

%Outputs: 
%histwgt: nbinsx1 vector: share of weights in each bin 
%midbin: nbinsx1: midpoint value for each bin

delta = (endpt-startpt)/nbins; 
subs = ceil((vv-startpt)/delta); 
subs(subs==0)=1;
midbin=((startpt+ delta/2):delta:endpt)';
histwgt = accumarray(subs,ww,[nbins,1]); 
end


 


 
