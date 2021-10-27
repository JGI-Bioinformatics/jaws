/* eslint-disable react-hooks/exhaustive-deps */
import React, { useRef, useEffect, useState } from 'react';

import * as d3 from "d3";

import pipeline2workflow from "../../../utils/pipeline_parser";
import classes from "./d3styles.module.css";
import randomId from "../../../utils/strings";

const NodeClassMap = {
  completed : classes.completed,
  failed : classes.failed,
  running : classes.running,
}

const Dag = (props) => {
  // console.log(props, "props, [Dag] render")
  const {width, height, data, charttype} = props;
  // const [svgdata, setSVGData] = useState();
  const [nodes, setNodes] = useState([]);
  // const [edges, setEdges] = useState([])

  const ref = useRef();

  const css = useRef({
    dragline : classes.dragline,
    hidden : classes.hidden,
    link : classes.link
  })

  const state = useRef({
    selectedG: null,          // selected "g" object containing "circle"

    selectedEdge: null,
    justDragged: false,
    justScaleTransGraph: false,
    lastKeyDown: -1,
    shiftNodeDrag: false,
    selectedText: null,
    pressedKey: null
  });

  const consts = useRef({
    selectedClass: "selected",
    connectClass: "connect-node",
    circleGClass: "conceptG",
    graphClass: "graph",
    activeEditId: "active-editing",
    radius: 40,
    width : width ? width : 400,
    height : height ? height : 400,
    BACKSPACE_KEY: 32,
    DELETE_KEY: 8,
    ENTER_KEY: 13,
    SHIFT_KEY: 16,
  }) 
 
  const THIS = useRef({
    nodes: [],
    edges: [],
    idct: 0,
    svgG: null,
    paths: null,
    circles: null
  });

  // display or hide tooltip
  const dispalyTip = (d=null, hide=true) => {
    // console.log(d, hide, "[displayTip]")
    const ttip = d3.select(ref.current).select(".node-tooltip");
    if (hide) {
      ttip.style("display", "none");
    } else if (d && d.title){
      ttip.style("left", d3.event.pageX - 30 + "px")
          .style("top", d3.event.pageY - 30 + "px")
          .style("display", "inline-block")
          .html(d.title);
    }
  }

  const resetGraph = () => {
    // console.log(THIS.current.nodes, THIS.current.edges, "[resetGraph]")
    let paths = getPaths();
    
    paths = paths.data(THIS.current.edges, d=>d);
      // update existing edges
      paths.style('marker-end','url(#end-arrow)')
      .classed(classes.link, true)
      .attr("d", d => `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`);

      // remove old links
      paths.exit().remove();
    
    // update new edlges
    paths.enter()
          .append("path")
          .style('marker-end','url(#end-arrow)')
          .classed(classes.link, true)
          .attr("d", d => `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`);

    // conditionally install event handler on edges
    if (charttype === "build"){
      getPaths().data(THIS.current.edges, d=>`${d.source.id}->${d.target.id}`)
                .on("dblclick", function(d){
                  d3.event.stopPropagation();
                  // console.log("dbclick on edge");
                  setSelectedEdge(d3.select(this));
                })
                .on("mousedown", function(d){
                  // console.log("mousedown on edge");
                  setSelectedEdge(d3.select(this));
                })
                .on("mouseover", function(d){
                  // console.log(d3.select(this), "mouseover on edge")
                })
    }
          

    const circles = getCircles().data(THIS.current.nodes, d => d.id);

    // add new nodes;
    const nodeSVG = circles.enter()
      .append("g")
      .attr("transform", d =>  `translate(${d.x},${d.y})`)
      .attr("id", d => d.id)
      .call(THIS.current.drag);
    
    // conditionally install event handler on edges
    if (charttype === "build"){
      getCircles().data(THIS.current.nodes, d => d.id)
                  .on("dblclick", function(d){
                    // console.log("dbclick on node");
                    d3.event.stopPropagation();
                    doubleClickOnNode(d3.select(this));
                  })
    } else if (charttype === "abstract"){
      getCircles().data(THIS.current.nodes, d => d.id)
                .on("mouseover", (d)=>{
                  // console.log("mouseover on node");
                  d3.event.stopPropagation();
                  dispalyTip(d, false);
                })
                .on("mouseout",() => { 
                  // console.log("mouseout on node");
                  d3.event.stopPropagation();
                  dispalyTip();
                })
                
    }
      

      
    // now nodeSVG holds th newly added "g" list of the above!
    
    nodeSVG.append("circle")  // add circle to each of the newly added "g" of the above
      .attr("r", consts.current.radius)
      .attr("class", d => { return charttype !== "build" ? nodeDefaulClass(d) : classes.default})

    if (charttype === "build"){
      nodeSVG.each(function(d){     // add text to each of the newly added "g" of the above
        insertTitleLinebreaks(d3.select(this), d);
      });
    };

    // remove old nodes
    circles.exit().remove();
  };

  const renderGraph = () => {
    const paths = getPaths().data(THIS.current.edges);
    paths.style('marker-end', 'url(#end-arrow)')
      .attr("d", d => `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`);  

    getCircles()
        .data(THIS.current.nodes, d => d.id)
        .attr("transform", d => `translate(${d.x},${d.y})` );
  }

  const getPaths = () => (
    d3.select(ref.current).select("svg").select(".pathroot").selectAll("path")
  );

  const getCircles = () => (
    d3.select(ref.current).select("svg").select(".circleroot").selectAll("g")
  );

  // current mouse position in SVG, relative to SVG origin {x:INT, y:INT}
  const curMousePositionOnSVG = () => {
    return {x: d3.event.offsetX, y: d3.event.offsetY};
  };

  const nodeDefaulClass = data => {
    let cls = classes.idle;
    if (data && data.data && data.data.status){
      cls = NodeClassMap[data.data.status];
    }
    return cls;
  };


  /**
   * g : the new "g"
   */
  const setSelectedG = (g, d=null) => {
    deselectSelectedG();
 
    g.classed(classes.selected, true);
    state.current.selectedG = g;
  };

  // const setSelectedEdge = (e, d=null) => {

  //   deselectSelectedEdge();
  //   e.classed(classes.selected, true);
  //   state.current.selectedEdge = e;
  // };

  const deselectSelectedG = () => {
    const { selectedG } = state.current;
    if (selectedG){
      selectedG.classed(classes.selected, false);
      state.current.selectedG = null;
      THIS.current.dragLine.classed(css.current.hidden, true);
    };
  };

  const deselectSelectedEdge = () => {
    const { selectedEdge } = state.current;
    if (selectedEdge){
      selectedEdge.classed(classes.selected, false);
      state.current.selectedEdge = null;
    };
  };

  const doubleClickOnNode = (g) => {
    // console.log("[doubleClickOnNode]");
    const { selectedG } = state.current;
    const gdata = g.data()[0]
    if (selectedG){
      const sdata = selectedG.data()[0];
      if (sdata.id !== gdata.id) { // double click on a node different from previously dblClicked node
        deselectSelectedG();
        // console.log(" >>> add edge")
        addNewEdge(sdata, gdata);
      } else {
        deselectSelectedG();
      }
    } else {
      setSelectedG(g);
    }
  };

  const setSelectedEdge = (e) => {
    // console.log(e, e.data()[0], "[setSelectedEdge]");
    const { selectedEdge } = state.current;
    const edata = e.data()[0];
    let sameEdge = false;
    if (selectedEdge){
      const sdata = selectedEdge.data()[0];
      sameEdge = (edata.source.id === sdata.source.id && edata.target.id === sdata.target.id)
    };

    deselectSelectedEdge();
    
    if (!sameEdge) {
      e.classed(classes.selected, true);
      state.current.selectedEdge = e;
    };
  };

  const resetSelections = () => {
    THIS.current.dragLine.classed(css.current.hidden, true);
    state.current.selectedEdge = null;
    deselectSelectedG();
  };

  const addNewEdge = (fromNodeData, toNodeData) => {
  
    if (fromNodeData && toNodeData && fromNodeData !== toNodeData){
      // we're in a different node: create new edge for mousedown edge and add to graph
      const newEdge = {source: fromNodeData, target: toNodeData};
    
      // console.log(paths, "[addnewEdge] paths")
      const filtRes = THIS.current.edges.filter(function(d){
        if (d.source === newEdge.target && d.target === newEdge.source){
          THIS.current.edges.splice(THIS.current.edges.indexOf(d), 1);
        }
        return d.source === newEdge.source && d.target === newEdge.target;
      });

      if (filtRes.length === 0){
        // console.log("[addNewEdge ok]")
        THIS.current.edges.push(newEdge);
        resetGraph();
        renderGraph();
      }
    } 
      
    resetSelections();
  };

  /**
   * delete the given node object form THIS.current.nodes, together with any edges connecting to the node.
   */
  const deleteNode = (node) => {
    // remove the node from node list
    THIS.current.nodes.splice(THIS.current.nodes.indexOf(node), 1);

    // find all edges to the node
    const edges = THIS.current.edges.filter( e => 
                                            (e.source.id === node.id || e.target.id === node.id));
    // console.log(node, edges, "[deleteNode] to be rmoved")
    // remove all the connected edges
    edges.forEach( e => {
      THIS.current.edges.splice(THIS.current.edges.indexOf(e), 1);
    });

    resetGraph();
    renderGraph();
  };

  const deleteAction = () => {
    // console.log("[deleteAction]");
    const { selectedG, selectedEdge } = state.current;
    if ( selectedG ) {
      const ndata = selectedG.data()[0];
      const nodelist = THIS.current.nodes.filter(cd=>cd.id === ndata.id); // dobj 
      deleteNode(nodelist[0]);
      deselectSelectedG();       
    }

    if ( selectedEdge ) {
      const edata = selectedEdge.data()[0];
      deselectSelectedEdge();
      const edgelist = THIS.current.edges.filter(cd=>cd.source.id === edata.source.id && cd.target.id === edata.target.id);
      THIS.current.edges.splice(THIS.current.edges.indexOf(edgelist[0]), 1);
      
      resetGraph();
      renderGraph();
    }
  };

  // add a new node to data, update display and return the new node data;
  const addNewNode =() => {
    const { x, y } = curMousePositionOnSVG();
    const d = {id: THIS.current.idct++, title: "new step", x, y};
    THIS.current.nodes.push(d);
    resetGraph();
   
    // make title of text immediently editable
    const d3txt = changeTextOfNode(getCircles().filter( dval => dval.id === d.id), d);
    const txtNode = d3txt.node();
    selectElementContents(txtNode);
    txtNode.focus();   
    return d;
  };

  const dragmove = d => {
    // console.log(state.current, "dragmove")
    d3.select(ref.current).select(".node-tooltip").style("display", "none");
    if (state.current.shiftNodeDrag){
      // console.log(`M${d.x},${d.y}L${d3.mouse(THIS.current.svgG.node())[0]},${d3.mouse(THIS.current.svgG.node())[1]}`, "[mousemove]")
      THIS.current.dragLine.attr('d', `M${d.x},${d.y}L${d3.mouse(THIS.current.svgG.node())[0]},${d3.mouse(THIS.current.svgG.node())[1]}`);
    } else {
      d.x += d3.event.dx;
      d.y += d3.event.dy;
      renderGraph();
      //console.log(`(${d3.event.dx}; ${d3.event.dy}) event`)
    }
  };

  const insertTitleLinebreaks = (gEl, d) => {
    // console.log(gEl, d, "insertTitleLinebreaks")
    const words = d.title.split(/\s+/g);
    const nwords = words.length;
    let yoffset = (nwords-1)*7.5;
    // if (nwords === 1) { yoffset -= 40 }

    const el = gEl.append("text")
          .attr("text-anchor","middle")
          .attr("dy", "-" + yoffset)
          .classed(classes.Blue, true)
          .on("dblclick", ()=>{         //double click on label to edit
            // console.log("doubleClickOnText ...")
            d3.event.stopPropagation();
            const d3txt = changeTextOfNode(gEl, d);
            const txtNode = d3txt.node();
            selectElementContents(txtNode);
            txtNode.focus();
          })

    for (let i = 0; i < nwords; i++) {
      const tspan = el.append('tspan').text(words[i]);
      if (i > 0)
        tspan.attr('x', 0).attr('dy', '15');
    }
  };

  /* select all text in element: taken from http://stackoverflow.com/questions/6139107/programatically-select-text-in-a-contenteditable-html-element */
  const selectElementContents = el => {
    const range = document.createRange();
    range.selectNodeContents(el);
    const sel = window.getSelection();
    sel.removeAllRanges();
    sel.addRange(range);
  };

  /* place editable text on node in place of svg text */
  const changeTextOfNode = function(d3node, d){
    const htmlEl = d3node.node();
    
    // console.log(d, "changeTextOfNode")
    d3node.selectAll("text").remove();
    const nodeBCR = htmlEl.getBoundingClientRect(),
        curScale = nodeBCR.width/consts.current.radius,
        useHW = curScale > 1 ? nodeBCR.width*0.71 : consts.current.nodeRadius*1.42;

    // replace with editableconent text
    // var d3txt = THIS.current.svg.selectAll("foreignObject")
    const d3txt = d3node.selectAll("foreignObject")
          .data([d])
          .enter()
          .append("foreignObject")
          .attr("x", -consts.current.radius )
          .attr("y", -consts.current.radius + 12)
          .attr("height", 2*useHW)
          .attr("width", 2*consts.current.radius)
          .append("xhtml:p")
          .attr("id", consts.activeEditId)
          .attr("contentEditable", "true")
          .text(d.title)
          .classed(classes.editP, true)
          .on("mousedown", d => { d3.event.stopPropagation();})
          .on("keydown", function(d) {
            d3.event.stopPropagation();
            if (d3.event.keyCode === consts.current.ENTER_KEY && !d3.event.shiftKey){
              this.blur();
            }
          })
          .on("blur", function(d){
            d.title = this.textContent;
            insertTitleLinebreaks(d3node, d);
            d3.select(this.parentElement).remove();
          });
          
    return d3txt;
  };

  const svgKeyDown = function() {
    console.log("[svgKeyDown]");
    if(state.current.lastKeyDown !== -1) return;
    state.current.lastKeyDown = d3.event.keyCode;
  };

  const svgKeyUp = () => { 
    // console.log(d3.select("svg").select(".pathroot").selectAll("path"), "[svgKeyUp]");
    console.log(d3.event.keyCode, "[svgKeyUp]")
    state.current.lastKeyDown = -1;

    if(d3.event.keyCode === consts.current.DELETE_KEY){
      deleteAction();
    }

    /* seems when keydown, the mouseUp event didn't get trigger
     */
    if (state.current.shiftNodeDrag){
      // dragged from node
      state.current.shiftNodeDrag = false;
      THIS.current.dragLine.classed(css.current.hidden, true);
    }
  };


  useEffect( () => {
    // console.log("[Dag] useEffect default")
    const radius_arrow_map = {
      50: 80,
      40: 66,
      25: 45,
      15: 32,
      10: 22
    };    
    
    consts.current.radius = charttype === "abstract" ? 10 : consts.current.radius;
    
    const svg = d3.select(ref.current).append('svg')
      .attr("width", consts.current.width)
      .attr("height", consts.current.height)
      .style("border", "1px solid black");

    // define arrow markers for graph links
    
    const refx = `${radius_arrow_map[consts.current.radius]}`;
    THIS.current.svg = svg;
    svg.append('defs')
      .append('marker')
      .attr('id', 'end-arrow')
      .attr('viewBox', '0 -5 10 10')
      .attr('refX', refx)
      .attr('markerWidth', 3.5)
      .attr('markerHeight', 3.5)
      .attr('orient', 'auto')
      .append('path')
      .attr('d', 'M0,-5L10,0L0,5')
      // define arrow markers for leading arrow (drag between 2 circles)
      .append('marker')
      .attr('id', 'mark-end-arrow')
      .attr('viewBox', '0 -5 10 10')
      .attr('refX', refx)
      .attr('markerWidth', 3.5)
      .attr('markerHeight', 3.5)
      .attr('orient', 'auto')
      .append('path')
      .attr('d', 'M0,-5L10,0L0,5');
    
    const topclass = randomId(12); // need this be unique in case of multiple instances in one page: zoom and translate
    const svgG = svg.append("g").classed(consts.current.graphClass, true).classed(topclass, true);

    // displayed when dragging between nodes
    THIS.current.dragLine = svgG.append('path')
      // .attr('class', 'link dragline hidden')
      .classed(css.current.link, true)
      .classed(css.current.dragline, true)
      .classed(css.current.hidden, true)
      .attr('d', 'M0,0L0,0')
      .style('marker-end', 'url(#mark-end-arrow)');

    svgG.append("g").classed("pathroot", true).selectAll("g");
    svgG.append("g").classed("circleroot", true).selectAll("g");
    THIS.current.svgG = svgG;

    THIS.current.drag = d3.drag()
      .on("start", d => ({x: d.x, y: d.y}))
      .on("drag", args => {
        THIS.current.justDragged = true;
        dragmove(args);
      })


    if (charttype === "build"){
      d3.select(window)
      .on("keydown", ()=>{ svgKeyDown();})
      .on("keyup", () => { svgKeyUp();});

      svg.on("mousemove", () => {
        // console.log(d3.event, "[mousemove]")
        if (state.current.selectedG) {  // draw the dragLine
          const { x, y } = curMousePositionOnSVG();
          const d = state.current.selectedG.select("circle").data()[0];
          const { dragLine } = THIS.current;
          dragLine.classed(css.current.hidden, false);
          dragLine.attr('d', `M${d.x},${d.y}L${x},${y}`);
        }
      })
      .on("dblclick", function(){
        //d3.event.stopPropagation();
        // console.log("svgDoubleClick");
        if(state.current.lastKeyDown !== consts.current.SHIFT_KEY){
          const newNode = addNewNode();
          
          // if has selectedG, create edge
          const { selectedG } = state.current;
          if ( selectedG !== null ) {
            const sdata = state.current.selectedG.select("circle").data()[0];
            addNewEdge(sdata, newNode);
            // deselectSelectedG();
          }
        }
      })
    }
    
    if (charttype === "abstract" ){ // zoom() interfers with mouseup event !
      // listen for dragging
      const dragSvg = d3.zoom()
      .scaleExtent([0.8, 4])   // zoom range
      .on("zoom", () => {
        if (d3.event.sourceEvent.shiftKey){
          // TODO  the internal d3 state is still changing
          return false;
        } else{
          
          // justScaleTransGraph = true;
          const currentTransform = d3.event.transform;
          d3.select(`.${topclass}`)
            .attr("transform", currentTransform); 
        }
        return true;
      })
      svg.call(dragSvg).on("dblclick.zoom", null);
    }

    // window
    // d3.select(window)
    // .on("keydown", function(){
    //   thisGraph.svgKeyDown.call(thisGraph);
    // })
    // .on("keyup", function(){
    //   thisGraph.svgKeyUp.call(thisGraph);
    // });
    // svg.on("mousedown", function(d){thisGraph.svgMouseDown.call(thisGraph, d);});
    // svg.on("mouseup", function(d){thisGraph.svgMouseUp.call(thisGraph, d);});


    if (charttype === "build" && process.env.NODE_ENV !== "production"){
      // debug:
      d3.select("#debug").selectAll("button")
        .data(["nodes", "edges", "resetGraph", "reRender", "state", "toggleDragLine"])
        .enter()
        .append("button")
        .text(d=>d)
        .on("click", d => {
          if (d === "nodes"){
            // console.log(THIS.current.nodes, `[DEBUG] ${d}`);
          } else if (d === "edges") {
            // console.log(THIS.current.edges, `[DEBUG] ${d}`);
          } else if (d === "resetGraph") {
            resetGraph();
            // console.log(`[DEBUG] ${d}`);
          } else if (d === "reRender"){
            // console.log(`[DEBUG] ${d}`);
          } else if (d === "state"){
            // console.log(state.current, `[DEBUG] ${d}`)
          } else if (d === "toggleDragLine"){
            THIS.current.dragLine.classed(css.current.hidden, !THIS.current.dragLine.classed(css.current.hidden));
          }
        })
    }
  }, []);

  useEffect(() => {
    // console.log("[Dag] useEffect props")

    let xLoc = 60;
    let xSpace = 120;
    let ySpace = 100;
    if (charttype === "abstract"){
      xLoc = 15;
      xSpace = 46;
      ySpace = 24;
      // consts.current.radius = 10;
    }
    
    const yLoc = consts.current.height / 2;
    const wfdata = pipeline2workflow(data, xLoc, yLoc, xSpace, ySpace);

    // console.log(wfdata, "set cvdata")
    THIS.current.idct = wfdata.nodes.length;
    THIS.current.nodes = wfdata.nodes;
    THIS.current.edges = wfdata.edges;  
    // console.log(wfdata, "[Dag] useEffect on data")
    // setSVGData(wfdata);
    setNodes(wfdata.nodes);
    // setEdges(wfdata.edges);
 
  }, [data, charttype]);

  useEffect(() =>{
    // console.log(nodes, "[dag] useEffect on wfdata")
    if (nodes){
      resetGraph();
    }
  }, [nodes]);

  
  const style = [];
  if ( charttype === "abstract" ){
    style.push(classes.ToolTip);
    style.push("node-tooltip");
  };

  return (
    <>
      <div ref={ref}>
        <div className={style.join(" ")}></div>
      </div>
      <div id="debug">

      </div>
    </>
  )
}

export default Dag;
