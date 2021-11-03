
import { withRouter } from 'react-router';
// import { withStyles } from '@material-ui/core/styles';
// import Table from '@material-ui/core/Table';
// import TableBody from '@material-ui/core/TableBody';
// import TableCell from '@material-ui/core/TableCell';
// import TableContainer from '@material-ui/core/TableContainer';
// import TableHead from '@material-ui/core/TableHead';
// import TableRow from '@material-ui/core/TableRow';
// import Paper from '@material-ui/core/Paper';

// import classes from "./MyJobs.module.css";


// const StyledTableCell = withStyles((theme) => ({
//     head: {
//       backgroundColor: theme.palette.common.black,
//       color: theme.palette.common.white,
//     },
//     body: {
//       fontSize: 14,
//     },
//   }))(TableCell);
  
//   const StyledTableRow = withStyles((theme) => ({
//     root: {
//       '&:nth-of-type(odd)': {
//         backgroundColor: theme.palette.action.hover,
//       },
//     },
//   }))(TableRow);

function Run(props) {
    console.log(props, "[Run]");
    let record = JSON.parse(localStorage.getItem('fetchedData'))[props.match.params.i];
    /*const more_info = 
    <TableContainer className={classes.tableBox} component={Paper}>
      <Table className={classes.table} aria-label="simple table">
        <TableHead>
        <TableRow>
            <StyledTableCell></StyledTableCell>
            <StyledTableCell>More Information on {record["id"]}</StyledTableCell>
        </TableRow>
        </TableHead>
        <TableBody>
        <TableRow>
            <StyledTableCell>Status Detail</StyledTableCell>
            <StyledTableCell>{record["status_detail"]}</StyledTableCell>
        </TableRow>
        <TableRow>
            <StyledTableCell>WDL</StyledTableCell>
            <StyledTableCell>{record["wdl_file"]}</StyledTableCell>
        </TableRow>
        <TableRow>
           <StyledTableCell>JSON</StyledTableCell>
            <StyledTableCell>{record["json_file"]}</StyledTableCell>
        </TableRow>
        <TableRow>
            <StyledTableCell>Tag</StyledTableCell>
            <StyledTableCell>{record["tag"]}</StyledTableCell>
        </TableRow>
        </TableBody>
        </Table>
    </TableContainer>*/
    
    const header = ( 
        record["status_detail"] === "" &&
        record["wdl_file"] == null &&
        record["json_file"] == null &&
        record["tag"] === null) ?
        <h4>
            No more information on {record["id"]}!
        </h4>:
        <h4>
            More information on {record["id"]}:<br/>
        </h4>
    
    const status_detail = record["status_detail"] !== "" ? 
        <>Status Detail: {record["status_detail"]}</> : <></>
    
    const wdl_file = record["wdl_file"] != null ? 
        <>WDL File: {record["wdl_file"]}</> : <></>
    
    const json_file = record["json_file"] != null ? 
        <>JSON File: {record["json_file"]}</> : <></>
    
    const tag = record["tag"] != null ? 
        <>Tag: {record["tag"]}</> : <></>
    
    
    return (
        <div>
            {header}
            {status_detail}
            {wdl_file}
            {json_file}
            {tag}
        
        </div>
    )

    /*
    return (
        <div>
            {more_info}
        </div>
    )*/
}

export default withRouter(Run);