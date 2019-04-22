function[fillhandle,msg]=jbfill(xpoints,upper,lower,roi,color,edge,add,transparency) 
%% ADAPTED BY K. OLFERS 28-08-2016
% https://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves
%USAGE: [fillhandle,msg]=jbfill(xpoints,upper,lower,color,edge,add,transparency) 
%This function will fill a region with a color between the two vectors provided 
%using the Matlab fill command. 
% 
%fillhandle is the returned handle to the filled region in the plot. 
%xpoints= The horizontal data points (ie frequencies). Note length(Upper) 
% must equal Length(lower)and must equal length(xpoints)! 
%upper = the upper curve values (data can be less than lower) 
%lower = the lower curve values (data can be more than upper) 
%roi = a *(logical) vector specifying the areas to be shaded 
%color = the color of the filled area 
%edge = the color around the edge of the filled area 
%add = a flag to add to the current plot or make a new one. 
%transparency is a value ranging from 1 for opaque to 0 for invisible for 
%the filled color only. 
% 
%John A. Bockstege November 2006; 
%Example: 
% a=rand(1,20);%Vector of random data 
% b=a+2*rand(1,20);%2nd vector of data points; 
% x=1:20;%horizontal vector 
% [ph,msg]=jbfill(x,a,b,rand(1,3),rand(1,3),0,rand(1,1)) 
% grid on 
% legend('Datr') 
if nargin<8;transparency=.4;end %default is to have a transparency of .5 
if nargin<7;add=0;end %default is to add to current plot 
if nargin<6;edge='k';end %dfault edge color is black 
if nargin<5;color='b';end %default color is blue 
if nargin<4; roi = ones(length(xpoints),1); end 
if ~isrow(lower); lower = lower'; end 
if ~isrow(upper); upper = upper'; end 
if ~isrow(xpoints); xpoints = xpoints'; end 
[clusters cnum] = bwlabeln(roi);
if ~add 
%figure; 
%plot([1:length(upper)],upper,[1:length(upper)],lower); 
add = 1; 
end 
for iClust = 1:cnum
iStart = find(clusters==iClust,1); 
clength = length(find(clusters == iClust))-1; 
tmpUpper = upper(iStart:iStart+clength); 
tmpLower = lower(iStart:iStart+clength); 
tmpXpoints = xpoints (iStart:iStart+clength); 

if length(tmpUpper)==length(tmpLower) && length(tmpLower)==length(tmpXpoints) 
msg=''; 
filled=[tmpUpper,fliplr(tmpLower)]; 
tmpXpoints=[tmpXpoints,fliplr(tmpXpoints)]; 
if add 
hold on 
end 
fillhandle=fill(tmpXpoints,filled,color);%plot the data 
ax = gca;
%ax.XTickLabel = xpoints;
datetick('x','dd-mmm-yyyy')
%h = gca; 
%ax.XAxis.TickLabelFormat = 'dd-mmm-yyyy';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color 
if add 
hold off 
end 
else 
msg='Error: Must use the same number of points in each vector'; 
end
end