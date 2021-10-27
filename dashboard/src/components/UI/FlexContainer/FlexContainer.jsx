import React from 'react';
import CssBaseline from '@material-ui/core/CssBaseline';
import classes from "./FlexContainer.module.css";
import { Button } from "@material-ui/core";

import {cookieGet, CK_AVATAR_ID} from "../../../utils/cookie";

export default function FlexContainer(props) {
  let curId = cookieGet(CK_AVATAR_ID);
  curId = curId ? parseInt(curId) : 0;
  const [index, setIndex] = React.useState(curId);

  
  const items = props.items.map( (itm, idx) => { 
        let itmClass = [classes.ButtonLow]
        if (idx === index){
          // itmClass.push(classes.ButtonHigh);
          itmClass = [classes.ButtonHigh]
        }
        // console.log(style, idx, index, "[FlexContainer]")
        return ( <Button className={itmClass.join(" ")}
                          key={idx} 
                          onClick={ () => {setIndex(idx); props.click(idx)} }> 
                    {itm} 
                </Button>
                )
      })
  
  return (
    <React.Fragment>
      <CssBaseline />
      <div className={classes.FlexContainer}>
        {items}
      </div>
    </React.Fragment>
  );
}

