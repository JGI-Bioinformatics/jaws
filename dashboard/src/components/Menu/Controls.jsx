import React from 'react';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';

import { withRouter } from "react-router-dom";

import classes from "./Controls.module.css";

const VerticalTabs = (props) => {
  // console.log(props, "VerticalTabs")
  // const {flags, user} = props;
  // console.log(props, "[Control]")

  const tabDef = [...props.controls];
 
  const tabs = tabDef.map( (pair) => {
    const css = { display: pair[2]}
    return <Tab key={pair[0]} value={pair[1]} label={pair[0]} style={css}/>
  });

  // console.log(props, "DEBUG");
  return (
      <Tabs
        orientation="vertical"
        variant="scrollable"
        // textColor="primary"
        value={props.location.pathname}
        onChange={(e, value)=>{props.history.push(value)}}
        aria-label="Vertical tabs example"
        className={classes.Tabs}
      >
        {tabs}
      </Tabs>
  );
  
}

export default withRouter(VerticalTabs);