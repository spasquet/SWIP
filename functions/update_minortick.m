function update_minortick(minor_tick_no)

ax = gca;
ax.YMinorTick = 'on';
ax.YAxis.MinorTickValuesMode = 'manual'; % prevents MATLAB form update it
tick_gap = ax.YAxis.TickValues(2)-ax.YAxis.TickValues(1);

minor_gap = tick_gap/minor_tick_no;
ax.YAxis.MinorTickValues = ax.YAxis.TickValues(1)+minor_gap:...
    minor_gap:ax.YAxis.TickValues(end);