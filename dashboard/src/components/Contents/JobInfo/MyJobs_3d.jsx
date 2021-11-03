import React, {useState} from 'react';

import { connect } from "react-redux";

import BottomNavigation from '@material-ui/core/BottomNavigation';
import BottomNavigationAction from '@material-ui/core/BottomNavigationAction';
import DirectionsRunIcon from '@material-ui/icons/DirectionsRun';
import EmojiPeopleIcon from '@material-ui/icons/EmojiPeople';

// import FlowChart from "../../UI/Chart/FlowChart";
// import MyJobList from "./MyJobList";

// import pipe_data_1 from "./pipeline_data_1.js";
// import pipe_data_2 from "./pipeline_data_2.js";
// import user_jobs from "./myjob_data.js";

import MyJobActive from "./MyJobActive";
import MyJobPassed from "./MyJobOld";

import classes from "./MyJobs.module.css";



const MyJobs = props => {

  // use props.user.user_token to query user jobs
  const [content, setContentValue] = useState("active");

  const contentComp = content === "active" ? <MyJobActive /> : <MyJobPassed />

  const handleContentChange = (event, newValue) => {
    setContentValue(newValue);
  };

  return (
    <>
      {/* <FlowChart d3id="d3flowchart" width="1400" height="300" fromjob data={pipe_data_2}/> */}
      {/* <MyJobList page={page} /> */}
      {/* {props.user.email} */}
      <BottomNavigation value={content} onChange={handleContentChange} className={classes.Menu}>
        <BottomNavigationAction label="Active Jobs" value="active" icon={<DirectionsRunIcon />} />
        <BottomNavigationAction label="Finished Jobs" value="passed" icon={<EmojiPeopleIcon />} />
      </BottomNavigation>
      <div className={classes.Page}>
        { contentComp }
        
      </div>
    </>
  )
}

const mapStateToProps = state => {
  return {
    user: state.auth.user,
  }
}

export default connect(mapStateToProps)(MyJobs);