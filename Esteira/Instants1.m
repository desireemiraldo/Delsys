function [HeelStrike, ToeOff] = Instants1(instant)
% Select instants Heel Strike and Toe Off for the chosen trials
% index - indexes of the chosen trials


HeelStrike = instant(:,1:2:end);

ToeOff = instant(:,2:2:end);

end