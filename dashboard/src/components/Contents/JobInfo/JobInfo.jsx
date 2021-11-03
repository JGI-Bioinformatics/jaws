import React, { useState } from 'react';

import Divider from '@material-ui/core/Divider';

import Pagination from "../../UI/Pagination/Pagination";


// import classes from "./JobInfo.module.css";
const JobInfo = props => {
  const [page, setPage] = useState(1);  // default paginatiom to page 1

  // const handleChangePage = (event, newPage) => {
  const handleChangePage = (newPage) => {
    // console.log(page, newPage, "[MyJobActive] page")
    setPage(newPage);
  };

  return (
    <>
      
      <h3>Job information : page {page}</h3>
      <Divider />
      <Pagination initpage="1" count="10" onchange={handleChangePage} />
    </>
  )
}
 
export default JobInfo;