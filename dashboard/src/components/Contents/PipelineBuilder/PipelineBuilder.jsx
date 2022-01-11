import React, {useState, useEffect} from 'react';

// import { useState } from 'react';
// import FlowChart from "../../UI/Chart/FlowChart";
// import withLD from "../../hoc/withLD";

import Dag from "../../UI/Chart/Dag";
import data from "./data"

import classes from "./PipelineBuilder.module.css";

const PipelineBuilder = props => {
  const [wfdata, setWFData] = useState();

  useEffect(()=>{
    // console.log(props.flags.workflow, "[Workflow] useEffect once:wfdata in props");
    
    // if (props.flags && props.flags.workflow && Object.keys(props.flags.workflow).length > 0){
    //   setWFData(props.flags.workflow);  // wfdata from LD flag
    // } else {
      setWFData(data);
    // }
  }, [props])


  // useEffect(()=>{
  //   // log every rerender
  //   console.log("[Workflow] render")
  // })

  const content = wfdata ? <Dag d3id="d3flowchart" charttype="build" width={750} height={400} data={wfdata}/> :
                          <span> Loading ... </span>

  // const content = wfdata ? <FlowChart d3id="d3flowchart" width={600} height={400} data={wfdata}/> :
  //                         <span> Loading ... </span>

  const style = [classes.lightColor, classes.footNote];

  return ( 
    <div className="MainPanel">
      <h3> PipelineBuilder : </h3>
      {content}
      <hr />
      <div className={style.join(" ")}>
        <ul>
          <li>Select/Deselect Node/Edge : Double Click on object</li>
          <li>Select/Deselect Edge : SHIFT + MouseOver + Click</li>
          <li>Create Node : Double Click on a blank area (but NO SHIFT key down).</li>
          <li>Create Edge : Select a node and then Double Click on another node or a blank area.</li>
          <li>Delete Node/Edge : Press DELETE key to delete currently selected node / edge.</li>
        </ul>
      </div>
      <hr />
    </div>
  );
  
}
 
export default PipelineBuilder;