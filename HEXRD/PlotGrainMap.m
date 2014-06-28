function [] = PlotGrainMap(glist_log)

nGrains = length(glist_log);

figure,
for i = 1:1:nGrains
    plot3(glist_log(i).COM(1), glist_log(i).COM(2), glist_log(i).COM(3), 'ko')
    hold on
end