function mean_shaded_std(m_vector, std_vector)

x = 1:numel(m_vector);
% std_dev=1;
c1 = m_vector+std_vector;
c2 = m_vector-std_vector;

x2 = [x, fliplr(x)];
inBetween = [c1, fliplr(c2)];
h=fill(x2,inBetween,'g', 'edgecolor', 'none');
set(h, 'facealpha', .9)
alpha(.5)
hold on;
plot(x,m_vector,'r','LineWidth',2);
end