cbticklabsnow = compose('%g',linspace(cmin,cmax,5))';
if cmax<max(max(mapvarnow))
    cbticklabsnow{end} = ['>' cbticklabsnow{end}];
end
if cmin>min(min(mapvarnow))
    cbticklabsnow{1} = ['<' cbticklabsnow{1}];
end
cb.Ticks = linspace(cmin,cmax,5);
cb.TickLabels = cbticklabsnow;
