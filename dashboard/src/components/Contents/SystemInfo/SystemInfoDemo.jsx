import React, { Component } from "react";

import {Grid} from "@material-ui/core"
import BarChart from "../../UI/Chart/BarChart";
import DataTable from "../../UI/Table/DataTable";

import ClusterStatus from "./ClusterStatus/ClusterStatus";

import data from "./data";


class SystemInfo extends Component {
  state = {  }

  workerTableData = (data) => {
    //const header = [];
    const rows = [];
    for (const itm of data){
      // const cfg = itm.config;
      rows.push([
              {value: itm.name, align: "left"},
              // {value: itm.total, align: "center"},
              {value: itm.busy, align: "center"},
              {value: itm.total - itm.busy, align: "center"},
            ]);
    }
    return rows;
  }

  configTableData = (data) => {
    const rows = [];
    for (const itm of data){
      const cfg = itm.config;
      rows.push([
              {value: itm.name, align: "left"},
              {value: itm.total, align: "center"},
              {value: `${cfg.memory.value}${cfg.memory.unit}`, align: "center"},
              {value: `${cfg.timelimit.value}${cfg.timelimit.unit}`, align: "center"},
              {value: cfg.cpu, align: "center"},
            ]);
    }
    return rows;
  }

  render() { 
    const configRows = this.configTableData(data.workers);
    const configHeader = [
              {label: "Worker Type", align: "left"},
              {label: "Count", align: "center"},
              {label: "Memory", align: "center"},
              {label: "Time Limit", align: "center"},
              {label: "# CPU", align: "center"},
            ];
    // const workerRows = this.workerTableData(data.workers);
    // const workerHeader = [
    //             {label: "Type", align: "left"},
    //             // {label: "Total", align: "center"},
    //             {label: "Busy", align: "center"},
    //             {label: "Idle", align: "center"},
    //           ]
    return ( 
      <div className="MainPanel">
        <h3> System Configuration and Status. </h3>
        <Grid container spacing={2} justify="flex-start">
          <Grid item xs={12} sm={6} md={8}  >
            <BarChart data={data.workers}></BarChart>
            {/* <DataTable rows={workerRows} header={workerHeader} title="Worker Status"/> */}
          </Grid>
          <Grid item xs={12} sm={4} md={4} >
            <div style={{width: "200px", height: "200px", paddingLeft: "20px"}}> 
              <ClusterStatus cluster={data.jtm} />
            </div>
          </Grid>
        </Grid>
        
        <DataTable rows={configRows} header={configHeader} title="Configuration"/>
        {/* <DataTable rows={workerRows} header={workerHeader} title="Worker Status"/> */}
       
      </div>
    );
  }
}
 
export default SystemInfo;