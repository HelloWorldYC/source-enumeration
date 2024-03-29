% 自定义配色3，31种推荐的精致配色

clear;
clc;
close all;

rgbTriplet = [0.57, 0.69, 0.30;...
    0.89, 0.88, 0.57;...
    0.76, 0.49, 0.58;...
    0.47, 0.76, 0.81;...
    0.21, 0.21, 0.35;...
    0.28, 0.57, 0.54;...
    0.07, 0.35, 0.40;...
    0.41, 0.20, 0.42;...
    0.60, 0.24, 0.18;...
    0.76, 0.84, 0.65;...
    0.77, 0.18, 0.78;...
    0.21, 0.33, 0.64;...
    0.88, 0.17, 0.56;...
    0.20, 0.69, 0.28;...
    0.26, 0.15, 0.47;...
    0.83, 0.27, 0.44;...
    0.87, 0.85, 0.42;...
    0.85, 0.51, 0.87;...
    0.99, 0.62, 0.76;...
    0.52, 0.43, 0.87;...
    0.00, 0.68, 0.92;...
    0.26, 0.45, 0.77;...
    0.98, 0.75, 0.00;...
    0.72, 0.81, 0.76;...
    0.77, 0.18, 0.78;...
    0.28, 0.39, 0.44;...
    0.22, 0.26, 0.24;...
    0.64, 0.52, 0.64;...
    0.87, 0.73, 0.78;...
    0.94, 0.89, 0.85;...
    0.85, 0.84, 0.86];
len = size(rgbTriplet, 1);

figure(1);
L = 0.04; W = 0.92; B = 0.02; H = 0.94;
ax = axes('position', [L, B, W, H]);
set(gcf, 'unit', 'centimeters', 'position', [2 2 5 16])
x = 1:8;
y = log2(x);

axes(ax)
plot(x, y+1,...
    'color', rgbTriplet(1, :), 'LineWidth', 5);
hold on;
str = [handle_rgb(rgbTriplet(1, 1)),...
    '-', handle_rgb(rgbTriplet(1, 2)),...
    '-', handle_rgb(rgbTriplet(1, 3))];
text(8, log2(8)+1, str,...
    'FontName', 'Times New Roman',...
    'FontSize', 9)
for icolor = 2:len
    plot(x, y+icolor,...
        'color', rgbTriplet(icolor, :), 'LineWidth', 5);
    str = [handle_rgb(rgbTriplet(icolor, 1)),...
        '-', handle_rgb(rgbTriplet(icolor, 2)),...
        '-', handle_rgb(rgbTriplet(icolor, 3))];
    text(8, log2(8)+icolor, str,...
        'FontName', 'Times New Roman',...
        'FontSize', 9)
end
hold off;
ax.XLim = [1, 12];
ax.YLim = [1, len+log2(8)];
title('28 Kinds of Delicate Color');
ax.Box ='on';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
set(ax, 'FontSize', 10.5, 'FontName', 'Times New Roman');

print(gcf, '31_discrete_colors_recommended_in_matlab.png', '-dpng', '-r600')

function str = handle_rgb(rgbTriplet)
str = num2str(rgbTriplet);
str = str(2:end);
if isempty(str)
    str = '.00';
else
    if length(str) == 2
        str = [str, '0'];
    end
end

end