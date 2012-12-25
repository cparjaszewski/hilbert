% PLOTAX.M
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.
%
% PLOTAX.M	- Adds a (black) axis-cross to the existing plot. 
%
% Input		: color
%
% Output	: a axis-cross to the existing plot
%
% Usage		: plotax(color);
%    		  plotax;
%
% See also	: PLOTC, PLOTNUM, PLOTRI, PLOTRIS, PLOTST

% HGFei     : Nov. 1999
% Modified  : Jun. 2004, Vlad Vicol

function plotax(color)

    if nargin ==0; color = 'k'; end;   

    ax = axis; 

    if(length(ax)==4)
    hold; 
         plot(linspace(ax(1),ax(2),100), zeros(1,100),color);
         plot(zeros(1,100), linspace(ax(3),ax(4),100),color); 
    hold off;  

    end
    if(length(ax)==6)
    hold;
         plot3(linspace(ax(1),ax(2),100), zeros(1,100),zeros(1,100),color);
         plot3(zeros(1,100), linspace(ax(3),ax(4),100),zeros(1,100),color); 
         plot3(zeros(1,100),zeros(1,100), linspace(ax(5),ax(6),100),color); 
    hold off;
    end
    if(length(ax)==2)
    hold;
        plot(linsace(ax(1),ax(2),100),color);
    hold off;
    end


    figure(gcf);
end