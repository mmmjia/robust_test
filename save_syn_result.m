scale_all(1).s=1:3;
scale_all(2).s=4:8;
scale_all(3).s=9:11;
scale_all(4).s=13:14;
scale_all(5).s=15:18;
scale_all(6).s=19:22;
scale_all(7).s=24:43;

bt=[12,23];

timeline=8:50;

for i=1:size(scale_all,2)
scale=scale_all(i).s;
x=timeline(scale);
y=maxf(scale);
err=(varx1(scale)-vaex2(scale))/2;
hold on
fill([x,flip(x)],[varx1(scale),flip(vaex2(scale))],[0.322,0.847,0.859], 'EdgeColor', 'none', 'FaceAlpha', 0.5)
end

plot(timeline,maxf,'Color',[0.65 0.85 0.9, 1])


for i=1:size(scale_all,2)
scale=scale_all(i).s;
x=timeline(scale);
y=maxf(scale);
err=(varx1(scale)-vaex2(scale))/2;
hold on
errorbar(x,y,err,"-s","MarkerSize",10,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],"MarkerEdgeColor",[0.65 0.85 0.90])
end

for i=1:size(scale_all,2)
figure
timeline=scale_all(i).s+7;
for tt=timeline
fill(-0.5*save_voi0(tt-7).f0+tt,-save_voi0(tt-7).x,[0,1,1])
hold on
end
xlabel('Month','Fontsize',15)
ylabel('Frequency','Fontsize',15)

%{
ticky=yticks;
yticklabels(-ticky)

xticks(timeline)
tick=month(timeline-7);
xticklabels(tick)
box off
saveas(gcf, ['my_figure',num2str(i),'.jpg']);
close()
%}
end

for i=1:22
x=x_all(i).x;
f=f_all(i).f;
fill(-0.5*f+i,x,[0,1,1])
hold on
end
%%save single amino acid plot
