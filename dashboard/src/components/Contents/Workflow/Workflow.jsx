import React from "react";
import ExpandableTable from "./expandableTable";
import PageContainer from "../PageContainer";

const Workflow = (props) => {
 
  return ( 
    <PageContainer title="Workflow Catalog:" content={<ExpandableTable />}/>
  );
  
}
 
export default Workflow;