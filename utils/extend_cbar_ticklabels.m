cbticklabsnow = compose('%g',linspace(cmin,cmax,cbnumticks))';
if cmax<max(max(mapvarnow))
    cbticklabsnow{end} = ['>' cbticklabsnow{end}];
end
if cmin>min(min(mapvarnow))
    cbticklabsnow{1} = ['<' cbticklabsnow{1}];
end
cb.Ticks = linspace(cmin,cmax,cbnumticks);
cb.TickLabels = cbticklabsnow;
