function [newAxis] = LinkAxisData (XorY,TickPositions,TickLabels,AxisName)
% [newAxis] = LinkAxisData (XorY,TickPositions,TickLabels,AxisName)
%
% XorY specifies which axis to add to ('x' or 'y')
% 
% TickPositions specifies the position of the ticks for the top/right x/y axis in
% the units of the bottom/left axis.
%
% TickLabels specifies the numerical label for these ticks
%
% AxisName is optional, but if specified it will label the new axis
%
% newAxis is a handle to the newly created axis
%
% Example (Arrhenius plot with temperature as a top scale):
% T=300:-10:100;E=exp(1./T); % Define a dummy set of data
% figure;plot(1000./T,E);xlabel('1000/T / K^{-1}') % Plot the data
% LinkTopAxisData(1000./(300:-50:100),300:-50:100,'T / K'); % Add a top axis
%
% The function will create the new axis behind the current one, so
% clicking gives the original axis rather than the new (dummy) one.
% The properties of the new axis are also linked to the old one, so zooming
% and rescaling should work as advertised.  It will even handle changes to
% linear / log etc.
%
% The (known) bugs are:
%
%  1) On a y log scale where one of the limits is zero the tick marks won't
%  quite align - this is because Matlab reports the lower limit as zero
%  even though it displays a non zero lower limit
%
%  2) Because the two axes are linked, the autoscale checkbox doesn't work
%  for the x and y limits.
%
%  This code is based on the addTopXAxis code by Emmanuel P. Dinnat
%  which can be found on the Matlab file exchange.  The main differences
%  between his program and this one are
%
%  1) This code inputs the top ticks as an array, and they don't have to
%  correspond to the tick positions on the bottom axis
%
%  2) This code hides the new axis behind the first one (making clicking
%  more intuitive)
%
%  3) However... this code is less powerful.
%
% Modifief from LinkTopAxisData.m, 20140528 GMW

XorY = lower(XorY);

oldAxis=gca;
set(get(oldAxis,'Title'),'Units','Normalized','Position',[0.5 1.05 0]); % Shift up the title a bit
% set(oldAxis,'Color','none','box','off');

% Create new axis
newAxis = axes('position', get(oldAxis, 'position'));
set(newAxis,'Color',get(oldAxis,'Color'));
set(newAxis, 'xGrid', 'off', 'yGrid', 'off'); % remove grids
set(newAxis,'XDir',get(oldAxis,'XDir'),'YDir',get(oldAxis,'YDir'))
    
switch XorY
    case 'x'
        % put new Xaxis on top
        set(newAxis, 'xaxisLocation', 'top');
        set(oldAxis, 'xaxislocation', 'bottom');

        % Give the new axis a name if necessary
        if nargin==4; xlabel(AxisName); end
        
        % Get rid of y labels on the new axis
        set(newAxis, 'yTickLabel', []);
        % simulate the presence of a box by making the new axis' y axis appear
        % on the right whilst the old one sits on the left
        set(newAxis, 'YAxisLocation', 'right');
        set(oldAxis, 'YAxisLocation', 'left');
        
        set(newAxis, 'YScale', get(oldAxis, 'YScale'));
        set(newAxis, 'YTick', get(oldAxis, 'YTick'));
        set(newAxis, 'XScale', get(oldAxis, 'XScale'));
        set(newAxis, 'XTick', get(oldAxis, 'XTick'));
        
        linkaxes([oldAxis newAxis],'xy');
        % Create a link between properties
        hlink = linkprop([newAxis, oldAxis], {'Position','YTick','XScale','YScale','YMinorTick','XDir','YDir'});
        % And store it on the new axis (to make sure it gets updated, but is
        % also destroyed when the figure is closed / axis is deleted)
        setappdata(newAxis,'Axis_Linkage',hlink);
        
        % Label the new axis bits
        set(newAxis,'XTick',TickPositions)
        set(newAxis,'XTickLabel',TickLabels)
        
    case 'y'
         % put new Xaxis on top
        set(newAxis, 'yaxisLocation', 'right');
        set(oldAxis, 'yaxislocation', 'left');
        % Give the new axis a name if necessary
        if nargin==4;ylabel(AxisName);end;
        
        % Get rid of y labels on the new axis
        set(newAxis, 'xTickLabel', []);

        % simulate the presence of a box by making the new axis' y axis appear
        % on the right whilst the old one sits on the left
        set(newAxis, 'XAxisLocation', 'top');
        set(oldAxis, 'XAxisLocation', 'bottom');
        
        set(newAxis, 'YScale', get(oldAxis, 'YScale'));
        set(newAxis, 'YTick', get(oldAxis, 'YTick'));
        set(newAxis, 'XScale', get(oldAxis, 'XScale'));
        set(newAxis, 'XTick', get(oldAxis, 'XTick'));
        
        linkaxes([oldAxis newAxis],'xy');
        % Create a link between properties
        hlink = linkprop([newAxis, oldAxis], {'Position','XTick','XScale','YScale','XMinorTick','XDir','YDir'});
        % And store it on the new axis (to make sure it gets updated, but is
        % also destroyed when the figure is closed / axis is deleted)
        setappdata(newAxis,'Axis_Linkage',hlink);
        
        % Label the new axis bits
        set(newAxis,'YTick',TickPositions)
        set(newAxis,'YTickLabel',TickLabels)
end
        
% and finally, swap the places of the two axes so that clicking gives the correct
% behaviour
temp=get(gcf,'Children');
i=temp==newAxis;
j=temp==oldAxis;
temp(i)=oldAxis;
temp(j)=newAxis;

set(gcf,'Children',temp);

axis(oldAxis);
