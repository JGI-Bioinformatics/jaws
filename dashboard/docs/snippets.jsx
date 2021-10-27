// import React from "react";
import React, {useEffect} from "react";

import React, { Component } from 'react';

class BlastDB extends Component {
  state = {  }
  
  componentDidMount() {

  }
  render() { 
    return (  );
  }
}
 
export default BlastDB;
const SystemInfo = () => {
  
  useEffect(() => {

      // httpclient.get(JAWS_API_STATUS).then(resp => {
      //   // const userinfo = resp.data;
      //   console.log(resp.data);
      // });

      const tid = setTimeout(()=>{
        console.log('I am done');
      }, 2000);

      return () => { clearTimeout(tid) }
  },[])

  console.log('[SystemInfo]')
  return <h1>This is my header</h1>
  //http://jaws.lbl.gov:5002/api/v2/status
  /*
  {
  "CORI-Cromwell": "UP",
  "CORI-RMQ": "UP",
  "CORI-Site": "UP",
  "JAWS-Central": "UP",
  "JGI-Cromwell": "UP",
  "JGI-RMQ": "UP",
  "JGI-Site": "UP"
}
*/


}

export default SystemInfo;