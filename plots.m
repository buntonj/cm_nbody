limits = [-0.5;1;-0.5;1;-0.5;1];
limits = (limits)*0.5E-8;
lastrun = importdata('lasttype.txt');	%Which potential data?
pottype = cell2mat(lastrun.textdata);
pnum = lastrun.data(1);
tsteps = lastrun.data(2);
dogif = 1;	%Do the gif?

if pottype == 'L'
fprintf('Plotting L-J Potential \n')
filename = 'LJposition.gif';
else
fprintf('Plotting Morse Potential \n')
filename = 'Mposition.gif';
end
R = importdata('pos.txt');
tsteps = 1000;

initial = figure;
scatter3(R(1:pnum,2),R(1:pnum,3),R(1:pnum,4));
title('Initial Orientation');
xlabel('x (meters)');
ylabel('y (meters)');
zlabel('z (meters)');
if pottype == 'L'
saveas(initial,'LJinitial','epsc');
else
saveas(initial,'Minitial','epsc');
end

final = figure;
scatter3(R((tsteps)*pnum + 1:(tsteps+1)*pnum,2),R(tsteps*pnum + 1:(tsteps+1)*pnum,3),R(tsteps*pnum + 1:(tsteps+1)*pnum,4));
title('Final Orientation');
xlabel('x (meters)');
ylabel('y (meters)');
zlabel('z (meters)');
if pottype == 'L'
saveas(final,'LJfinal','epsc');
else
saveas(final,'Mfinal','epsc');
end

if dogif == 1
h = figure;
for i = 1:tsteps
%figure()
scatter3(R(((i-1)*pnum + 1):(i*pnum),2),R(((i-1)*pnum + 1):(i*pnum),3),R(((i-1)*pnum + 1):(i*pnum),4));
title(num2str(i))
axis(limits)
drawnow
% Capture the plot as an image
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if i == 1
imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
else
imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
end
end
end

E = importdata('E.txt');
P = importdata('P.txt');
V = importdata('vel.txt');

eplot = figure;
plot(E(:,1),E(:,2));
title('Total Energy vs. Time');
xlabel('time (seconds)');
ylabel('Energy (Joules)');
if pottype == 'L'
saveas(eplot,'LJEvsT','epsc');
else
saveas(eplot,'MEvsT','epsc');
end

pplot = figure;
plot(P(:,1),P(:,5));
title('Total Momentum vs. Time');
xlabel('time (seconds)');
ylabel('Momentum (kg*m/s)');
if pottype == 'L'
saveas(pplot,'LJPvsT','epsc');
else
saveas(pplot,'MPvsT','epsc');
end

vhist = figure;
histogram(V(1:125,2),20);
title('Initial Velocity Distribution');
xlabel('Velocity (m/s)');
ylabel('Frequency (count)');
saveas(vhist,'initialvelocity','epsc');

rave = zeros(tsteps,2);
for i = 1:tsteps
rave(i,2) = sum(R(((i-1)*pnum + 1):(i*pnum),5))/pnum;
rave(i,1) = R((i-1)*pnum + 1,6);
end
raver = figure;
plot(rave(:,1),rave(:,2));
title('Average Displacement vs. Time');
xlabel('time (seconds)');
ylabel('Displacement (meters)');
if pottype == 'L'
saveas(raver,'LJaverager','epsc');
else
saveas(raver,'Maverager','epsc');
end

adj = zeros(pnum,pnum);
for i = 1:pnum
for j = 1:pnum
	adj(i,j) = sqrt((R(tsteps*pnum+i,2)-R(tsteps*pnum+j,2))^2 + (R(tsteps*pnum+i,3)-R(tsteps*pnum+j,3))^2 + (R(tsteps*pnum+i,4)-R(tsteps*pnum+j,4))^2);
end
end
figure()
histogram(adj);
title('Final Relative Displacement');
xlabel('r (meters)');
ylabel('Frequency');
