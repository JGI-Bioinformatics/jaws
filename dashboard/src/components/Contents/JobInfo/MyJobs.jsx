import React, { useState } from "react";
import { connect } from "react-redux";
import { useQuery } from "react-query";
import { LinearProgress } from '@material-ui/core';
import CheckCircleOutlineIcon from '@material-ui/icons/CheckCircleOutline';
import TablePagination from '@material-ui/core/TablePagination';
import { Route, Link } from "react-router-dom";

import classes from "./MyJobs.module.css";
import MyJobsControl from "./Controls/MyJobsControl";
import httpclient, { JAWS_API_SEARCH } from "../../../utils/HttpClientProvider";
import DataTable from "../../UI/Table/DataTable";
import PageContainer from "../PageContainer";
import Run from "./Run";

// import mockjobs from "./mockjobs.json";


const fetchMyJobs = async (active, delta, site, result, access_token) => {
  const fdata = new FormData();
  fdata.append("active_only", active ? "True" : "False");
  fdata.append("delta_days", delta);
  fdata.append("site_id", site);
  fdata.append("result", result);

  // console.log(active, delta, site, result, access_token, "[fetchMyJobs]");

  const api = JAWS_API_SEARCH;
  const headers = {Authorization: `Bearer ${access_token}`};
  const resp = await httpclient.post(api, fdata, {headers: headers})
  // console.log(active, delta, site, result, resp.data, "[fetchMyJobs]");
  return resp.data;
  // return mockjobs;
}

const DELTA = "delta";
const SITE = "site";
const ACTIVE = "active";
const RESULT = "result";

const getDefault = (key, default_value) => {
  var value = localStorage.getItem(key);
  
  // boolean value is stored in string format.
  return value === null ? default_value : key === ACTIVE ? value === "true" : value;
}

const MyJobs = ({access_token}) => {

  // TODO : move into MyJobsControl  
  const [active, setActive] = useState(getDefault(ACTIVE, false)); 
  //check with Ed sice backed is takig a strig value instead of a boolea
  const [delta, setDelta] = useState(getDefault(DELTA, 100));
  const [site, setSite] = useState(getDefault(SITE, "all"));
  const [result, setResult] = useState(getDefault(RESULT,  "any"));
  
  const {data, status, error} = useQuery(
    ["myjobs", active, delta, site, result],
    () => fetchMyJobs(active, delta, site, result, access_token)
  );
  const [page, setPage] = useState(1);
  const [rowsPerPage, setRowsPerPage] = React.useState(5);

  const handleChangePage = (event, newPage) => {
    setPage(newPage);
  };

  const handleChangeRowsPerPage = (event) => {
    setRowsPerPage(parseInt(event.target.value, 10));
    setPage(0);
  };
  
  // console.log(active, delta, site, result, "[MyJobs]");

  const setSearchParams = (active, delta, site, result) => {
    // console.log(active, delta, site, result, "[setSearchParams]");
    // TODO : add a doFetch state and call update here 
    setActive(active);
    localStorage.setItem(ACTIVE, active);
    setDelta(delta);
    localStorage.setItem(DELTA, delta);
    setSite(site);
    localStorage.setItem(SITE, site);
    setResult(result);
    localStorage.setItem(RESULT, result);
  }

  if (status !== "success"){
    return <h2> loading ...</h2>;
  }

  const records_meta = data ? 
    <div className={classes.infoBox}> 
      <span className={classes.infoText}>
        <u><b>{data.length}</b></u> records found</span>
      <span 
        id='checkSign' 
        style={{
          verticalAlign: 'middle',
          marginLeft: '5px'}}> 
        <CheckCircleOutlineIcon className={classes.infoCheckIcon} fontSize="small"/></span>
        
    </div> : 
    error ? `ERROR : ${error}` : <LinearProgress />
  
  const header = [
    {label: "ID"}, 
    {label: "Input Site"},
    {label: "Site"},
    {label: "Result"},
    {label: "Status"},
    {label: "Submitted"},
    {label: "Updated"},
    {label: "More"},
  ];

  const rows = data ? data
                          .slice(page * rowsPerPage, page * rowsPerPage + rowsPerPage)
                          .map((row, idx) => {
                            return [
                              {
                                value: row.id 
                              },
                              {
                                value: row.input_site_id, 
                              },
                              {
                                value: row.site_id, 
                              },
                              {
                                value: row.results, 
                              },
                              {
                                value: row.status, 
                              },
                              {
                                value: row.submitted, 
                              },
                              {
                                value: row.updated, 
                              },
                              {
                                value: <Link to={"/run/" + ((page * rowsPerPage) + idx)}>
                                read more</Link>,
                              },
                            ]
                          }) : [];
  
  // TODO : move into MyJobsControl itself 
  const contralData = {
                        delta, 
                        site,
                        result,
                        setAll: setSearchParams,
                        active,
                      }

  const content = <>
                    <MyJobsControl {...contralData}/>
                    {records_meta}
                    <Route exact path="/run/:i"><Run/></Route>
                    { data?.length > 0 ? <DataTable header={header} rows={rows}/>: null }
                    {
                      data?.length > 0 && 
                          <TablePagination
                            rowsPerPageOptions={[5, 10, 25]}
                            component="div"
                            count={data.length}
                            rowsPerPage={rowsPerPage}
                            page={page}
                            onChangePage={handleChangePage}
                            onChangeRowsPerPage={handleChangeRowsPerPage}
                          />
                    }
                  </>

  return (
    <PageContainer 
      title="My jobs :" 
      content={content} />
  );
}
const mapStateToProps = state => {
  return {
    access_token: state.auth.access_token,
  }
}
export default connect(mapStateToProps)(MyJobs);