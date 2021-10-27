import React from 'react';
import CircularProgress from '@material-ui/core/CircularProgress';
import { useQuery } from "react-query";
import ReadmeContent from "./readmeContent"

import DataTable from "../../UI/Table/DataTable";


// //fetch('https://code.jgi.doe.gov/qaqc/jgi-rqc-pipeline/-/blob/master/README.md')
// curl -X GET --header "Accept: application/json" --header "Authorization: Bearer LLKFALD..." "https://code.jgi.doe.gov/api/v4/groups/136"
const fetchWDLCatalog = async () => {
  const resp = await fetch("https://code.jgi.doe.gov/api/v4/groups/136", {  
     method: 'GET',
     headers: {
       'Authorization': `Bearer ${process.env.REACT_APP_MY_PWD}`,
     }
  });
  return await resp.json();
}

export default function ExpandableTable() {
  const { data, status } = useQuery("fetchWDLCatalog", fetchWDLCatalog)

  let content = null;

  const dateformatter = (d) => {
    // "2020-12-01T02:05:30.338Z" to "12/01/2020"
    const m = d.split("T").shift().split("-");
    return `${m[1]}/${m[2]}/${m[0]}`;
  }
  
  if (data?.message){
    // capture error case - can be { message: "401 Unauthorized" }
    content = data.message;
  } else if (status === "loading"){
    content = <CircularProgress />;
  } else {
    // Filter out project ids that I don't want included 
    // And add a new key "shared".
    const data_projects = data.projects.filter((row) => (![174,177].includes(row.id)) )

    // // adding key "shared" to "shared_projects"
    // data_projects.forEach(e => e['shared']='original project')
    // data.shared_projects.forEach(e => e['shared']='shared link')

    // Combine the two datasets from the shared and un-shared projects.
    const final_data = data_projects.concat(data.shared_projects);
    // console.log(final_data, "[expandableTable]")

    const header = [
                      {label: "Repository Name"}, 
                      {label: "Summary"},
                      {label: "Last Modified"},
                      {label: "Owner"},
                      {label: "Readme"},
                    ];
    const rows = final_data.map((row) => { 
        return [
          {
            value: <a href={row.web_url} target="_blank" rel="noreferrer">{row.name}</a> 
          },
          {
            value: row.description ? row.description : "none", 
          },
          {
            value: dateformatter(row.last_activity_at), 
          },
          {
            value: row.owner ? row.owner.name : "unknown", 
          },
          {
            value: <div>
                   <ReadmeContent repo_id={row.id} /> 
                   </div>
          },
        ]
      })

    content = <DataTable header={header} rows={rows}/>;
  }
  return (
    <div>
      {content}
    </div>
  );
}

