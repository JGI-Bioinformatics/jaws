import React from "react";
import Box from "@material-ui/core/Box";

import Search from "./Search/Search";
import SystemInfo from "./SystemInfo/SystemInfo";
import JobInfo from "./JobInfo/JobInfo"
import Workflow from "./Workflow/Workflow"
import User from "./User/User"

const ContentPanel  = (props) => {
  const pageMap = [
    <Search />,
    <SystemInfo />,
    <JobInfo/>,
    <Workflow/>
  ]

  if (props.useridentity){
    pageMap.unshift(<User useridentity={props.useridentity}/>)
  }

  return (
    <div>
        <Box p={3}>
          {pageMap[props.pageId]}
        </Box>
    </div>
  );
}

export default ContentPanel;