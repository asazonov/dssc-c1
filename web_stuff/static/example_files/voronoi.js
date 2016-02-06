function draw() {
  var bbox = d3.select(".jumbotron").node().getBoundingClientRect();

  var width = bbox.width,
      height = bbox.height;

  var vertices = d3.range(100).map(function(d) {
    return [Math.random() * width, Math.random() * height];
  });

  var voronoi = d3.geom.voronoi()
      .clipExtent([[0, 0], [width, height]]);

  d3.select(".voronoi").selectAll("*").remove();
  var svg = d3.select(".voronoi").append("svg")
      .attr("width", width)
      .attr("height", height)
      .on("mousemove", function() { vertices[0] = d3.mouse(this); redraw(); });

  var path = svg.append("g").selectAll("path");

  svg.selectAll("circle")
      .data(vertices.slice(1))
    .enter().append("circle")
      .attr("transform", function(d) { return "translate(" + d + ")"; })
      .attr("r", 0);

  redraw();

  function redraw() {

    path = path
        .data(voronoi(vertices), polygon);

    path.exit().remove();

    path.enter().append("path")
        .attr("class", function(d, i) { return "q" + (i % 9) + "-9"; })
        .attr("d", polygon);

    path.order();
  }

  function polygon(d) {
    return "M" + d.join("L") + "Z";
  }
}
draw();
window.addEventListener('resize', function(event){
    draw();
});

var myVar = setInterval(myTimer, 1000000);

function myTimer() {
    draw()
}
