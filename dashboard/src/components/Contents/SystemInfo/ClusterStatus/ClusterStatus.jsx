import React, { Component } from 'react';

import classes from "./ClusterStatus.module.css";

class ClusterStatus extends Component {
  state = {  }
  render() { 
    const rows = this.props.cluster.map( (itm, idx) => {
              const color = itm.up ? classes.CircleGreen : classes.CircleRed;                ;
              return (
                <tr key={idx}>
                  <td className={classes.Box}><div className={color}></div></td>
                  <td title={itm.description}>{itm.name}</td>
                </tr>
              )
            })
    return ( 
      <>
        <h4 className={classes.Header}>JTM Status</h4>
        <table>
          <tbody>
            {rows}
          </tbody>
        </table>
      </>
     );
  }
}
 
export default ClusterStatus;