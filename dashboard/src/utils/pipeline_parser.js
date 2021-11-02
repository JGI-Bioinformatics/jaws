/*
  take take pipeline json and return workflow chart ready nodes and edges;
  pdata = {
    steps: [
      { "name": STEP_NAME_STR,   // required
        "level": POSITION_IN_WORKFLOW_NUM(0-based), // required
        "branch": {       // parallel steps - optional
          "size": NUM_OF_PARAM_JOBS_IN_THIS_STEP_NUM, 
          "position": -1,  // INDEX AROUND 0 FOR size : size:2 -> -1, 1; size:3 -> -1,0,1; size:4 -> -2,-1,1,2 etc; for position node around level
        },
        "data": {
          "status": "completed" // STATUS_STR, completed | inprogress | fail ete; required
          "start": TIME_STR;
          "end": TIME_STR;
          //any other data
        }
      }
    ];
    workflow: [];
  }
*/

const pipeline2workflow = (pdata, xstart, ystart, xspace, yspace) => {
  const xLoc = xstart;
  const xSpace = xspace;
  const yLoc = ystart;

  const nodes = pdata.steps.map( (n, idx) => { 
    let y = yLoc;
    if (n.branch){
      const ySpace = yspace;
      y += n.branch.position * ySpace;
    }

    const stepobj =  {title: n.name, 
                        id: idx, 
                        x: xLoc + xSpace * n.level, 
                        y
                      };
   
    stepobj.data = n.data;
    
    return stepobj;
  })

  const edges = pdata.workflow.map( e => {
    const edge = {...e};

    pdata.steps.forEach( (n, idx) => {
      if (edge.source === n.name){
        edge.source = nodes[idx];
      }
      if (edge.target === n.name){
        edge.target = nodes[idx]
      }
    })
    return edge;
  });

  return {nodes, edges};
};

export default pipeline2workflow


