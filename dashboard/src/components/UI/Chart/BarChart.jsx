import React, { useEffect, useCallback, useRef } from "react";

import * as d3 from "d3"

import classes from "./BarChart.module.css";

/**
 * data : [{name:STR, busy:INT, total:INT}]
 * id: STR the dive 
 */
const BarChart = ({ data }) => {
 
  const ref = useRef();

  const drawChart = useCallback((data, subgroups) => {
    // console.log(data, "[Data]")
    // console.log(subgroups, "[subgroups]")

    const groups = data.map( d => d.name );
    // console.log(groups, "[groups]")

    const yRange = Math.max(...data.map(d=>d.idle+d.busy));
 
    // set the dimensions and margins of the graph
    const margin = {top: 40, right: 10, bottom: 40, left: 50},
        width = 460 - margin.left - margin.right,
        height = 400 - margin.top - margin.bottom;

    // append the svg object to the body of the page
    const svg = d3.select(ref.current)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);
    
    const tooltip = d3.select(ref.current).append("div").attr("class", classes.ToolTip);
    
    // Add X axis
    const x = d3.scaleBand()
                .domain(groups)
                .range([0, width])
                .padding([0.4]);

    svg.append("g")
      .attr("transform", `translate(0,${height})`)
      .style("font-size", 18)
      .call(d3.axisBottom(x).tickSizeOuter(0));

    // Add Y axis
    const y = d3.scaleLinear()
                .domain([0, yRange])
                .range([ height, 0 ]);
                
    svg.append("g")
        .call(d3.axisLeft(y));

    // color palette = one color per subgroup
    const color = d3.scaleOrdinal()
      .domain(subgroups)
      .range(["#F00", "#999"])

    //stack the data? --> stack per subgroup
    const stackedData = d3.stack()
      .keys(subgroups)
      (data)

    console.log(stackedData, "[BarChart] stackedData")
    
    // Show the bars
    svg.append("g")
    .selectAll("g")
    // Enter in the stack data = loop key per key = group per group
    .data(stackedData)
    .enter().append("g")
    .attr("fill", d => color(d.key))
    .selectAll("rect")
    // enter a second time = loop subgroup per subgroup to add all rectangles
    .data( d => d )
    .enter().append("rect")
    .attr("x", d => x(d.data.name) )
    .attr("y", d => y(d[1]) )
    .attr("height", d => y(d[0]) - y(d[1]) )
    .attr("width",x.bandwidth())
    .on("mousemove", d => {
        tooltip
          .style("left", `${d3.event.pageX - 30}px`)
          .style("top", `${d3.event.pageY - 30}px`)
          .style("display", "inline-block")
          .html(`${Object.keys(d.data).find(key => d.data[key] === d[1] - d[0])} ${d[1] - d[0]}`);
    })
    .on("mouseout",() => { tooltip.style("display", "none");} );

    // legend
    const legend = svg.append('g')
    .attr('class', 'legend')
    .attr('transform', 'translate(20, -40)')
    .style("font-style", "italic");

    legend.selectAll('rect')
        .data(subgroups)
        .enter()
        .append('rect')
        .attr('x', 0)
        .attr('y', function(d, i){
            return i * 18;
        })
        .attr('width', 14)
        .attr('height', 14)
        .attr('fill', function(d, i){
            return color(i);
        });

    legend.selectAll('text')
        .data(subgroups)
        .enter()
        .append('text')
        .text(function(d){
            return d;
        })
        .attr('x', 18)
        .attr('y', function(d, i){
            return i * 18;
        })
        .attr('text-anchor', 'start')
        .attr('alignment-baseline', 'hanging');
  }, [])

  useEffect(() => {
    const cdata = data.map( d => ({ name: d.name, busy: d.busy, idle: d.total-d.busy}) );
    const subgroups = ["busy", "idle"];
    drawChart(cdata, subgroups);
  }, [data, drawChart])
  
  return <div ref={ref} ></div>  
}
 
export default BarChart;