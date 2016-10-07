export default function plotScatter(chart, data, x, y) {
  const { data: vals } = data;

  const dotUpdate = chart.selectAll('.dot').data(vals, d => d);
  dotUpdate.exit().remove();
  const dotEnter = dotUpdate.enter()
    .append('circle')
      .attr('class', 'dot');
  dotUpdate.merge(dotEnter)
    .attr('r', 3)
    .attr('cx', d => x(d[0]))
    .attr('cy', d => y(d[1]));
}
