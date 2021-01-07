function element_plotter(theta,beta,center,a)

rect_size = a*[cos(theta),cos(theta)];
corner_left = center - rect_size.*0.5;
hold on
colorvalue = theta*0.7+0.5;%(-1)^(i+j+1)+0.5;
colorvalue = min(0.85,colorvalue);
cm = jet;
colorID = max(1,sum(colorvalue > [0:1/length(cm(:,1)):1]));
myColor = cm(colorID, :);
set(gcf,'Color',[1,1,1]);
% h = rectangle('Position',[corner_left(1),corner_left(2),rect_size(1),rect_size(2)],...
%     'linewidth',1.5,'Curvature',3*a,'FaceColor',myColor);
points = [corner_left;corner_left+[rect_size(1),0];corner_left+rect_size;corner_left+[0,rect_size(2)]];
h = patch(points(:,1),points(:,2),myColor);
rotate(h,[0 0 1],beta*180/pi,[corner_left+0.5.*rect_size,0]);

lego_l = 18e-3;
lego_w = 3e-3;

Point1 = [-lego_l,-lego_w]';
Point2 = [-lego_l,+lego_w]';
Point3 = [-lego_w,+lego_w]';
Point4 = [-lego_w,+lego_l]';
Point5 = [+lego_w,+lego_l]';
Point6 = [+lego_w,+lego_w]';
Point7 = [+lego_l,+lego_w]';
Point8 = [+lego_l,-lego_w]';
Point9 = [+lego_w,-lego_w]';
Point10 = [+lego_w,-lego_l]';
Point11 = [-lego_w,-lego_l]';
Point12 = [-lego_w,-lego_w]';

for i=1:2
    for j=1:2
        thet = theta*(-1)^(i+j) - beta;
        RoMatrix = [cos(thet),sin(thet); -sin(thet),cos(thet)];
        Robeta = [cos(-beta),sin(-beta); -sin(-beta),cos(-beta)];
        u = Robeta*[i-1.5,j-1.5]'.*(cos(theta))*a/2 + center';
        Point(i,j,1,:) = u + (RoMatrix*Point1);
        Point(i,j,2,:) = u + (RoMatrix*Point2);
        Point(i,j,3,:) = u + (RoMatrix*Point3);
        Point(i,j,4,:) = u + (RoMatrix*Point4);
        Point(i,j,5,:) = u + (RoMatrix*Point5);
        Point(i,j,6,:) = u + (RoMatrix*Point6);
        Point(i,j,7,:) = u + (RoMatrix*Point7);
        Point(i,j,8,:) = u + (RoMatrix*Point8);
        Point(i,j,9,:) = u + (RoMatrix*Point9);
        Point(i,j,10,:) = u + (RoMatrix*Point10);
        Point(i,j,11,:) = u + (RoMatrix*Point11);
        Point(i,j,12,:) = u + (RoMatrix*Point12);
        Point(i,j,13,:) = u; % center
    end
end

 
for i=1:2
    for j=1:2
        X=[Point(i,j,1,1),Point(i,j,2,1),Point(i,j,3,1),...
            Point(i,j,4,1),Point(i,j,5,1),Point(i,j,6,1),...
            Point(i,j,7,1),Point(i,j,8,1),Point(i,j,9,1),...
            Point(i,j,10,1),Point(i,j,11,1),Point(i,j,12,1)];
        Y=[Point(i,j,1,2),Point(i,j,2,2),Point(i,j,3,2),...
            Point(i,j,4,2),Point(i,j,5,2),Point(i,j,6,2),...
            Point(i,j,7,2),Point(i,j,8,2),Point(i,j,9,2),...
            Point(i,j,10,2),Point(i,j,11,2),Point(i,j,12,2)];
        fill(X,Y,0.96.*[1,1,1],'LineWidth',1);
        hold on
        
%         plot([Point(i,j,13,1),Point(i,j,14,1)],...
%             [Point(i,j,13,2),Point(i,j,14,2)],'-y',...
%             'linewidth', 4);
        
%         axis([0,n+1,0,m+1])
    end
end
axis equal
axis off
end