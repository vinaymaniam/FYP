clear;
load(sprintf('pyHeirarchy4096_'));
h1 = heirarchy; 
i1 = index;
load(sprintf('pyHeirarchy4096'));
h2 = heirarchy; 
i2 = index;
load(sprintf('pyMap4096cell96_'));
map1 = Map;
load(sprintf('pyMap4096cell96'));
map2 = Map;
clear heirarchy index Map Res

% for i = 1:4096
%     m1 = map1{i};    
%     fprintf('[Good](max,min,mean) = %.3f   %.3f   %.3f\n',...
%         max(m1(:)),min(m1(:)),mean(m1(:)));
%     m2 = map2{i};
%     fprintf('[Bad](max,min,mean) = %.3f   %.3f   %.3f\n',...
%         max(m2(:)),min(m2(:)),mean(m2(:)));
%     fprintf('\n--------------------------------------------\n')
% end
%% REMARK- m1 is 0 mean for all i
for i = 1:64
    for k = 1:193
        m1 = h1(i,:,k);    
        fprintf('[Good](mean,std,median) = %.3f   %.3f   %.3f\n',...
            1e10*mean(m1(:)),std(m1(:)),median(m1(:)));
        m2 = h2(i,:,k);
        fprintf('[Bad](mean,std,median) = %.3f   %.3f   %.3f\n',...
            1e10*mean(m2(:)),std(m2(:)),median(m2(:)));
        fprintf('\n--------------------------------------------\n')
    end
end
%% REMARK- m2 has a much smaller mean than m1 for all i


% figure;
% p1=[]; p2 = [];
% for i = 1:64
%     for j = 1:193
%         p1 = [p1, h1(i,1,j)'];
%         p2 = [p2, h2(i,1,j)'];
%     end
% end
% plot(p1,'r.'); 
% hold on;
% plot(p2,'b.');


