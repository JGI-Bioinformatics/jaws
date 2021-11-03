import React, { Component } from 'react';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

import classes from "./DataTable.module.css"

class DataTable extends Component {
  /**
   * require props.header and props.rows
   * props.header : [ {label: STR, align: left|right}]
   * props.rows : [ [{value:ANY, align: left|right}], ]
   */
  state = {  }
  render() {
    const header = this.props.header.map( (h, index) => {
      // const align = index === 0 ? "left" : "right"
      return <TableCell 
              key={index} 
              align={h.align} 
              style={{width: h.width ? `${h.width}px` : "auto"}}
              className={classes.Header}>
                {h.label}
              </TableCell>
    });

    const rows = this.props.rows;

    return ( 
      <TableContainer className={classes.Main} component={Paper}>
        <div className={classes.Title}>{this.props.title}</div>
        <Table className={classes.Table} size="small" aria-label="a dense table">
          <TableHead>
            <TableRow>
              {header}
            </TableRow>
          </TableHead>
          <TableBody>
            { rows.map( (row, idx) => (
                <TableRow key={idx}>
                  { row.map( (cell, idx) => (
                          idx === 0 ? <TableCell key={idx} component="th" scope="row" align={cell.align}>{cell.value}</TableCell> :
                          <TableCell key={idx} align={cell.align}>{cell.value}</TableCell>
                        ))
                  }
                </TableRow>
              )
            )}
          </TableBody>
        </Table>
      </TableContainer>
     
   );
  }
}
 
export default DataTable;
