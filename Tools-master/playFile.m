function obj = playFile(myfile);
   load(myfile);
   
   obj = audioplayer(y, Fs);
%    obj.TimerFcn = 'showSeconds';
%    obj.TimerPeriod = 1;
   
   play(obj);
end

% function showSeconds
%    disp('tick')
% end