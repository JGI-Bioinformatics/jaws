import React, { Component } from "react";

import * as d3 from "d3";

import pparser from "../../../utils/pipeline_parser";

import classes from "./d3styles.module.css";

// Ref : /Users/syao/MyDev/lang/web/d3/directG.html
/*
        <FlowChart 
          d3id= [ selector ID, must be unique if > 1 component is used ] 
          width="280" 
          height="80" 
          charttype="abstract" [ abstract: a scale-down instance with limited mouse interaction support ]
          data={} [the workflow object]
        />

        workflow data example :
        {
          "steps": [ 
            { "name": "Initialize", 
              "level": 0,
              "data": {"status": "completed"}
            },
            { "name": "Filter",
              "level": 1,
              "data": {"status": "completed"}
            },
            { "name": "Read QC",
              "level": 2,
              "branch": { "size": 3, "position": -1},
              "data": {"status": "completed"}
            },
            { "name": "Assembly",
              "level": 2,
              "branch": { "size": 3, "position": 0},
              "data": {
                "status": "completed",
              }
            },
            { "name": "This is a supperlongstep",
              "level": 2,
              "branch": { "size": 3, "position": 1},
              
              "data": {
              "status": "completed",
            }
            },
            {
              "name": "Post Analysis",
              "level": 3,
              "data": {"status": ""}
            },
            {
              "name": "Jamo Prep",
              "level": 4,
              "data": {
                "status": ""
              }
            }
          ],
          "workflow": [
              {"source": "Initialize", "target": "Filter"},
              {"source": "Filter", "target": "Assembly"},
              {"source": "Filter", "target": "Read QC"},
              {"source": "Filter", "target": "This is a supperlongstep"},
              {"source": "Assembly", "target": "Post Analysis"},
              {"source": "Read QC", "target": "Post Analysis"},
              {"source": "This is a supperlongstep", "target": "Post Analysis"},
              {"source": "Post Analysis", "target": "Jamo Prep"}
            ]
        }
*/
class FlowChart extends Component {
  status = {
    selectedNode: null,
    selectedEdge: null,
    mouseDownNode: null,
    mouseDownLink: null,
    justDragged: false,
    justScaleTransGraph: false,
    lastKeyDown: -1,
    shiftNodeDrag: false,
    selectedText: null
  };

  consts =  {
    selectedClass: "selected",
    connectClass: "connect-node",
    circleGClass: "conceptG",
    completed: "completed",
    failed: "failed",
    idle: "idle",
    graphClass: "graph",
    activeEditId: "active-editing",
    BACKSPACE_KEY: 8,
    DELETE_KEY: 46,
    ENTER_KEY: 13,
    nodeRadius: 40
  };

  //circle R to arrow offset map
  arrow_map = {
    50: 80,
    40: 66,
    25: 45,
    15: 32,
    10: 22
  }

 
  componentDidMount() {
    
    const {data, charttype} = this.props;
    // console.log(data, "[FlowChart]")

    // const charttype = "abstract";
    this.charttype = charttype;

    let xLoc = 60;
    let xSpace = 120;
    let ySpace = 100;
    if (charttype === "abstract"){
      xLoc = 15;
      xSpace = 46;
      ySpace = 24;
      this.consts.nodeRadius = 10;
    }
    
    this.setup(); // use the updated nodeRadius

    const yLoc = this.height / 2;
    const flowdata = pparser.pipeline2workflow(data, xLoc, yLoc, xSpace, ySpace);
    this.nodes = flowdata.nodes;
    this.edges = flowdata.edges;

    this.draw();    
  }

  setup = () => {
    const {d3id} = this.props;

    this.divid = `#${d3id}`;
    this.consts.graphClass = `${d3id}${this.consts.graphClass}`;
    this.canvas = d3.select(`#${d3id}`);
    // this.canvas.style("border", "1px solid #f00")

    let {width, height} = this.props;
    this.width = width ? width : 1200;
    this.height = height ? height : 600;
    // console.log(width, height)

    /** MAIN SVG **/
    this.svg = this.canvas.append("svg")
          .attr("width", this.width)
          .attr("height", this.height);

    this.tooltip = this.canvas.append("div").attr("class", classes.ToolTip);

    this.idct = 0;

    // define arrow markers for graph links
    const defs = this.svg.append('defs')
    const refx = `${this.arrow_map[this.consts.nodeRadius]}`;
    // console.log(refx, "[FlowChar] refx")
    defs.append('marker')
    .attr('id', 'end-arrow')
    .attr('viewBox', '0 -5 10 10')
    // .attr('refX', "32")
    .attr('refX', refx)
    .attr('markerWidth', 3.5)
    .attr('markerHeight', 3.5)
    .attr('orient', 'auto')
    .append('path')
    .attr('d', 'M0,-5L10,0L0,5')

    // define arrow markers for leading arrow (drag between 2 circles)
    defs.append('marker')
    .attr('id', 'mark-end-arrow')
    .attr('viewBox', '0 -5 10 10')
    .attr('refX', refx)
    .attr('markerWidth', 3.5)
    .attr('markerHeight', 3.5)
    .attr('orient', 'auto')
    .append('path')
    .attr('d', 'M0,-5L10,0L0,5');

    this.svgG = this.svg.append("g").classed(this.consts.graphClass, true).attr("id", "groot");
    
    // displayed when dragging between nodes
    this.dragLine = this.svgG.append('path')
      .attr('class', 'link dragline hidden')
      .attr('d', 'M0,0L0,0')
      .style('marker-end', 'url(#mark-end-arrow)');

    // svg nodes and edges 
    this.paths = this.svgG.append("g").classed("pathroot", true).selectAll("g");
    this.circles = this.svgG.append("g").classed("circleroot", true).selectAll("g"); 

    // listen for key events
    d3.select(window).on("keydown", ()=>{
      this.svgKeyDown();
    })
    .on("keyup", () => {
      this.svgKeyUp();
    });

    const svgMouseDown = this.svgMouseDown;
    const svgMouseUp = this.svgMouseUp;
    this.svg.on("mousedown", function(d){svgMouseDown(d);});
    this.svg.on("mouseup", function(d){svgMouseUp(d);});

    // if (this.charttype === "abstract"){
    //   return;
    // }

    // listen for dragging
    var dragSvg = d3.zoom()
          .scaleExtent([0.8, 4])   // zoom range
          .on("zoom", () => {
            if (d3.event.sourceEvent.shiftKey){
              // TODO  the internal d3 state is still changing
              return false;
            } else{
              this.zoomed();
            }
            return true;
          })
          // .on("zoomstart", () => {
          //   const ael = d3.select(`#${this.consts.activeEditId}`).node();
          //   if (ael){
          //     ael.blur();
          //   }
          //   if (!d3.event.sourceEvent.shiftKey) d3.select('body').style("cursor", "move");
          // })
          // .on("zoomend", () => {
          //   d3.select('body').style("cursor", "auto");
          // });
    
    this.svg.call(dragSvg).on("dblclick.zoom", null);

    // listen for resize
    window.onresize = ()=>{this.updateWindow(this.svg);};
  }

  updateGraph = () => {
    // console.log("[updateGraph]")
    d3.select(`${this.divid} .pathroot`).selectAll("path").data(this.edges)
      .attr("d", d => `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`);

    // update existing nodes
    d3.select(`${this.divid} .circleroot`).selectAll("g")
    .attr("transform", d => `translate(${d.x},${d.y})` );
  }

  draw(){
    // console.log("draw")

    const nodes = this.nodes;
    const edges = this.edges;
    const status = this.status;
    const consts = this.consts;

    // console.log(nodes, "[nodes]");
    // console.log(edges, "[edges");

    const paths = this.paths.data(edges, d => `${d.source.id}+${d.target.id}` ); // func don't have side effect!
    // const paths = this.paths.data(edges);

    // add new paths
    paths.enter()
      .append("path")
      .style('marker-end','url(#end-arrow)')
      .classed(classes.link, true)
      .attr("d", d => `M${d.source.x},${d.source.y}L${d.target.x},${d.target.y}`)
      // .on("mousedown", function(d){ pathMouseDown(d3.select(this), d); })
      // .on("mouseup", d => { status.mouseDownLink = null; });

    // remove old links
    paths.exit().remove();

    // update existing nodes
    this.circles = this.circles.data(nodes, d => d.id );
    this.circles.attr("transform", d => `translate(${d.x},${d.y})` );


    // add new nodes
    const circleMouseDown = this.circleMouseDown;
    const circleMouseUp = this.circleMouseUp;

    const newGs = this.circles.enter()
          .append("g")
          .classed(consts.circleGClass, d => !d.data)
          .attr("transform", function(d){return "translate(" + d.x + "," + d.y + ")";})
          .on("mouseover", function(d){        
            if (status.shiftNodeDrag){
              d3.select(this)
              .classed(consts.connectClass, true);
            }
          })
          .on("mouseout", function(d){
            d3.select(this).classed(consts.connectClass, false);
          })
          .on("mousedown", function(d){
            circleMouseDown(d3.select(this), d);
          })
          .on("mouseup", function(d){
            circleMouseUp(d3.select(this), d);
          })

          .on("mousemove", d => {
            if (d.title){
              this.tooltip
              .style("left", d3.event.pageX - 30 + "px")
              .style("top", d3.event.pageY - 30 + "px")
              .style("display", "inline-block")
              .html(d.title);
            }
          })
          .on("mouseout",() => { this.tooltip.style("display", "none");} )

          .call(d3.drag()
                .on("start", () => { 
                  this.tooltip.style("display", "none");
                  this.status.justDragged = true; 
                })
                .on("drag", args => { this.dragmove(args) })
                .on("end", () => { this.status.justDragged = false; }))

    newGs.append("circle")
      .attr("r", `${consts.nodeRadius}`)
      .classed(classes.completed, d => d.data && d.data.status === "completed")
      .classed(classes.failed, d => d.data && d.data.status === "failed")
      .classed(classes.idle, d => d.data && d.data.status === "")
      .classed(classes.running, d => d.data && d.data.status === "running")
 
    
    if (this.charttype !== "abstract"){
      const insertTitleLinebreaks = this.insertTitleLinebreaks
      newGs.each(function(d){
        insertTitleLinebreaks(d3.select(this), d);
      });
    };
  }

  dragmove = (d) => {
    // console.log(this.divid, "[dragmove]")
    if (this.status.shiftNodeDrag){
      this.dragLine.attr('d', 'M' + d.x + ',' + d.y + 'L' + d3.mouse(this.svgG.node())[0] + ',' + d3.mouse(this.svgG.node())[1]);
    } else{
      d.x += d3.event.dx;
      d.y +=  d3.event.dy;
      this.updateGraph();
      // console.log(`(${d.x}, ${d.y}) dragmove`)
    }
  };

  deleteGraph = skipPrompt => {
    let doDelete = true;
    if (!skipPrompt){
      doDelete = window.confirm("Press OK to delete this graph");
    }
    if(doDelete){
      this.nodes = [];
      this.edges = [];
      this.updateGraph();
    }
  };

  /* select all text in element: taken from http://stackoverflow.com/questions/6139107/programatically-select-text-in-a-contenteditable-html-element */
  selectElementContents = el => {
    const range = document.createRange();
    range.selectNodeContents(el);
    const sel = window.getSelection();
    sel.removeAllRanges();
    sel.addRange(range);
  };

  
  // remove edges associated with a node
  spliceLinksForNode = node => {
    const toSplice = this.edges.filter(function(l) {
      return (l.source === node || l.target === node);
    });

    toSplice.map( l => this.edges.splice(this.edges.indexOf(l), 1));
  };

  replaceSelectEdge = (d3Path, edgeData) => {
    d3Path.classed(this.consts.selectedClass, true);
    if (this.status.selectedEdge){
      this.removeSelectFromEdge();
    }
    this.status.selectedEdge = edgeData;
  };

  replaceSelectNode = (d3Node, nodeData) => {
    d3Node.classed(this.consts.selectedClass, true);
    if (this.status.selectedNode){
      this.removeSelectFromNode();
    }
    this.status.selectedNode = nodeData;
  };
  
  removeSelectFromNode = () => {
    this.circles.filter( cd => cd.id === this.status.selectedNode.id)
    .classed(this.consts.selectedClass, false);
    this.status.selectedNode = null;
  };

  removeSelectFromEdge = () => {
    this.paths.filter(cd => cd === this.status.selectedEdge)
    .classed(this.consts.selectedClass, false);
    this.status.selectedEdge = null;
  };

  pathMouseDown = (d3path, d) => {
    
    const status = this.status;
    d3.event.stopPropagation();
    status.mouseDownLink = d;

    if (status.selectedNode){
      this.removeSelectFromNode();
    }
    
    var prevEdge = status.selectedEdge;  
    if (!prevEdge || prevEdge !== d){
      this.replaceSelectEdge(d3path, d);
    } else{
      this.removeSelectFromEdge();
    }
  };

  insertTitleLinebreaks = (gEl, d) => {
    const words = d.title.split(/\s+/g);
    const nwords = words.length;
    let yoffset = (nwords-1)*7.5;
    // if (nwords === 1) { yoffset -= 40 }

    const el = gEl.append("text")
          .attr("text-anchor","middle")
          .attr("dy", "-" + yoffset)
          .classed("Blue", true)

    for (let i = 0; i < nwords; i++) {
      const tspan = el.append('tspan').text(words[i]);
      if (i > 0)
        tspan.attr('x', 0).attr('dy', '15');
    }
  };

  // mousedown on node
  circleMouseDown = (d3node, d) => {
    const status = this.status;
    d3.event.stopPropagation();
    status.mouseDownNode = d;
    if (d3.event.shiftKey){
      status.shiftNodeDrag = d3.event.shiftKey;
      // reposition dragged directed edge
      this.dragLine.classed('hidden', false)
        .attr('d', 'M' + d.x + ',' + d.y + 'L' + d.x + ',' + d.y);
      return;
    }
  };

  /* place editable text on node in place of svg text */
  changeTextOfNode = (d3node, d) => {
    const consts = this.consts;
    const insertTitleLinebreaks = this.insertTitleLinebreaks;

    const htmlEl = d3node.node();
    d3node.selectAll("text").remove();
    var nodeBCR = htmlEl.getBoundingClientRect(),
        curScale = nodeBCR.width/consts.nodeRadius,
        placePad  =  5*curScale,
        useHW = curScale > 1 ? nodeBCR.width*0.71 : consts.nodeRadius*1.42;
    // replace with editableconent text
    var d3txt = this.svg.selectAll("foreignObject")
          .data([d])
          .enter()
          .append("foreignObject")
          .attr("x", nodeBCR.left + placePad )
          .attr("y", nodeBCR.top + placePad)
          .attr("height", 2*useHW)
          .attr("width", useHW)
          .append("xhtml:p")
          .attr("id", consts.activeEditId)
          .attr("contentEditable", "true")
          .text(d.title)
          .on("mousedown", function(d){
            d3.event.stopPropagation();
          })
          .on("keydown", function(d){
            d3.event.stopPropagation();
            if (d3.event.keyCode === consts.ENTER_KEY && !d3.event.shiftKey){
              this.blur();
            }
          })
          .on("blur", function(d){
            d.title = this.textContent;
            insertTitleLinebreaks(d3node, d.title);
            d3.select(this.parentElement).remove();
          });
    return d3txt;
  };

  // mouseup on nodes
  circleMouseUp = (d3node, d) => {
    const status = this.status;
    const consts = this.consts;
    // reset the statuss
    status.shiftNodeDrag = false;    
    d3node.classed(consts.connectClass, false);
    
    const mouseDownNode = status.mouseDownNode;
    
    if (!mouseDownNode) return;

    this.dragLine.classed("hidden", true);

    if (mouseDownNode !== d){
      // we're in a different node: create new edge for mousedown edge and add to graph
      const newEdge = {source: mouseDownNode, target: d};
      const filtRes = this.paths.filter(function(d){
        if (d.source === newEdge.target && d.target === newEdge.source){
          this.edges.splice(this.edges.indexOf(d), 1);
        }
        return d.source === newEdge.source && d.target === newEdge.target;
      });
      if (!filtRes[0].length){
        this.edges.push(newEdge);
        this.updateGraph();
      }
    } else{
      // we're in the same node
      if (status.justDragged) {
        // dragged, not clicked
        status.justDragged = false;
      } else{
        // clicked, not dragged
        if (d3.event.shiftKey){
          // shift-clicked node: edit text content
          const d3txt = this.changeTextOfNode(d3node, d);
          const txtNode = d3txt.node();
          this.selectElementContents(txtNode);
          txtNode.focus();
        } else{
          if (status.selectedEdge){
            this.removeSelectFromEdge();
          }
          const prevNode = status.selectedNode;            
          
          if (!prevNode || prevNode.id !== d.id){
            this.replaceSelectNode(d3node, d);
          } else{
            this.removeSelectFromNode();
          }
        }
      }
    }
    status.mouseDownNode = null;
    return;
    
  }; // end of circles mouseup

  // mousedown on main svg
  svgMouseDown = () => {
    this.status.graphMouseDown = true;
  };

  // mouseup on main svg
  svgMouseUp = () => {

    const status = this.status;
    if (status.justScaleTransGraph) {
      // dragged not clicked
      status.justScaleTransGraph = false;
    } else if (status.graphMouseDown && d3.event.shiftKey){
      // clicked not dragged from svg
      const xycoords = d3.mouse(this.svgG.node()),
          d = {id: this.idct++, title: "new node", x: xycoords[0], y: xycoords[1]};
      this.nodes.push(d);
      this.updateGraph();
      // make title of text immediently editable
      const d3txt = this.changeTextOfNode(this.circles.filter(function(dval){
        return dval.id === d.id;
      }), d),
          txtNode = d3txt.node();
      this.selectElementContents(txtNode);
      txtNode.focus();
    } else if (status.shiftNodeDrag){
      // dragged from node
      status.shiftNodeDrag = false;
      this.dragLine.classed("hidden", true);
    }
    status.graphMouseDown = false;
  };

  // keydown on main svg
  svgKeyDown = () => {

    // const status = this.status;
    // const consts = this.consts;
    // // make sure repeated key presses don't register for each keydown
    // if(status.lastKeyDown !== -1) return;

    // status.lastKeyDown = d3.event.keyCode;
    // var selectedNode = status.selectedNode,
    //     selectedEdge = status.selectedEdge;

    // switch(d3.event.keyCode) {
    //   case consts.BACKSPACE_KEY:
    //   case consts.DELETE_KEY:
    //     d3.event.preventDefault();
    //     if (selectedNode){
    //       this.nodes.splice(this.nodes.indexOf(selectedNode), 1);
    //       this.spliceLinksForNode(selectedNode);
    //       status.selectedNode = null;
    //       this.updateGraph();
    //     } else if (selectedEdge){
    //       this.edges.splice(this.edges.indexOf(selectedEdge), 1);
    //       status.selectedEdge = null;
    //       this.updateGraph();
    //     }
    //     break;
    // }
  };

  svgKeyUp = () => {
    this.status.lastKeyDown = -1;
  };

  zoomed = () => {
    this.status.justScaleTransGraph = true;
    const currentTransform = d3.event.transform;
    d3.select("." + this.consts.graphClass)
      .attr("transform", currentTransform); 
  };

  updateWindow = svg => {
    const docEl = document.documentElement,
        bodyEl = document.getElementsByTagName('body')[0];
    const x = window.innerWidth || docEl.clientWidth || bodyEl.clientWidth;
    const y = window.innerHeight|| docEl.clientHeight|| bodyEl.clientHeight;
    svg.attr("width", x).attr("height", y);
  };



  render() { 
    const {d3id} = this.props;
    let {width, height } = this.props;
   
    return <div id={`${d3id}`} 
                style={ {width: `${width}px`, height: `${height}px`}} 
                className={classes.Chart}></div>  
  }
}
 
export default FlowChart;