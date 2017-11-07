function [Rf,Df] = splnwavf(wname) 
%     where the input argument wname is a string: 
%     wname = 'lem1' or 'lem2' ... i.e., 
%     wname = sh.name + number 
%     and w the corresponding scaling filter. 
%     The addition is obtained using:

if strcmp(wname,'spln4')
    Df=[1 4 6 4 1]./16;
    %Rf=[-5 20 -1 -96 70 280 70 -96 -1 20 -5]./256;%4
    %Rf=[3 -12 5 40 5 -12 3]./32;%3 to 4
    Rf=[35 -140 -55 920 -557 -2932 2625 8400 2625 -2932 -557 920 -55 -140 35]./8192; %5 to 4
    %Rf=[-63 189 329 -1568 -220 5948 -2732 -14752 10878 36750 10878 -14752 -2732 5948 -220 -1568 329 189 -63]./32768;
elseif strcmp(wname,'spln6')
    Df=[1 6 15 20 15 6 1]./64;
    Rf=[-63 378 -476 -1554 4404 1114 -13860 4158 28182 4158 -13860 1114 4404 -1554 -476 378 -63]./16384; %5 to 4
end
