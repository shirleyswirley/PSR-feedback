function [arrb,arrt] = computepointyColorbarlims(map,cmin,cmax)

pointiness=5; % smaller = pointier

if (min(map(:))<cmin)&(max(map(:))>cmax) % min and max are out of colorbar range
    arrb = cmin+(cmax-cmin)/pointiness;
    arrt = cmax-(cmax-cmin)/pointiness;
elseif min(map(:))<cmin % min is out of colorbar range
    arrb = cmin+(cmax-cmin)/pointiness;
    arrt = cmax;
    %arrt = cmax+(cmax-cmin); % start triangle even further off of colorbar
elseif max(map(:))>cmax % max is out of colorbar range
    arrb = cmin;
    %arrb = cmin-(cmax-cmin); % start triangle even further off of colorbar
    arrt = cmax-(cmax-cmin)/pointiness;
else % neither min nor max is out of colorbar range
    arrb = cmin; 
    arrt = cmax; 
    %arrb = cmin-(cmax-cmin); % start triangle even further off of colorbar
    %arrt = cmax+(cmax-cmin); % start triangle even further off of colorbar
end
