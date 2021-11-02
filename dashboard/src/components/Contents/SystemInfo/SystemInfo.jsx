import React from "react";

import { useQuery } from "react-query";
import CircularProgress from '@material-ui/core/CircularProgress';
import CheckCircleIcon from '@material-ui/icons/CheckCircle';
import CancelIcon from '@material-ui/icons/Cancel';

import { JAWS_API_STATUS } from "../../../utils/HttpClientProvider";
import DataTable from "../../UI/Table/DataTable";
import PageContainer from "../PageContainer";

import classes from "./SystemInfo.module.css";

const fetchSysInfo = async () => {
  const resp = await fetch(`https://jaws.lbl.gov:8011${JAWS_API_STATUS}`);
  return await resp.json();
}

const SystemInfo = ({jobid}) => {
  
  const { data, status } = useQuery("sysinfo",
                            fetchSysInfo
                          )

  let content = null;
  if (status !== "success"){
    content = <CircularProgress />;
  } else if (status === "success") {
    // console.log(data, "[data]");
    const header = [
      { label: "", width: "40"},
      { label: "Name"},
      { label: "Status", width: "100"},
    ];
    const rows = [];
    Object.keys(data).forEach( k => {
      rows.push([{value: data[k].toLowerCase() === "up" ? 
                          <CheckCircleIcon className={classes.upColor} /> : <CancelIcon className={classes.downColor} />,
                  align: "right" }, {value: k}, {value: data[k]}])
    })
    content = <DataTable header={header} rows={rows}/>;
  }

  return (
    <PageContainer 
    title="System Status :" 
    content={content} />
  )
}

export default SystemInfo;