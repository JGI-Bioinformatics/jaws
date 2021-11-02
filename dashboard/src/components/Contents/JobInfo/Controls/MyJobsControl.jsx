import React, { useState } from "react";
import { makeStyles } from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import MenuItem from '@material-ui/core/MenuItem';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Checkbox from '@material-ui/core/Checkbox';
import Button from '@material-ui/core/Button';
// import SearchIcon from '@material-ui/icons/Search';
import CheckBoxIcon from '@material-ui/icons/CheckBox';
// import Divider from '@material-ui/core/Divider';
import CheckBoxOutlineBlankIcon from '@material-ui/icons/CheckBoxOutlineBlank';
import styles from "./MyJobsControl.module.css";

const sites = [
  {value: 'cori', label: 'cori'},
  {value: 'jgi', label: 'jgi'},
  {value: 'all', label: 'all'},
];

const results = [
  {value: 'succeeded', label: 'succeeded'},
  {value: 'failed', label: 'failed'},
  {value: 'cancelled', label: 'cancelled'},
  {value: 'any', label: 'any'},
];

const useStyles = makeStyles((theme) => ({
  root: {
    '& .MuiTextField-root': {
      margin: theme.spacing(1),
    },
  },
}));

export default function JobsControl(props) {
  // console.log(props, "[MyJobsControl]props");
  const classes = useStyles();
  const [active, setActive] = useState(props.active);
  const [delta, setDelta] = useState(props.delta);
  const [site, setSite] = useState(props.site);
  const [result, setResult] = useState(props.result);

  const handleActive = (event) => {
    setActive(!active);
  };
  const handleDelta = (event) => {
    setDelta(event.target.value);
  };
  const handleSite = (event) => {
    setSite(event.target.value);
  };
  const handleResult = (event) => {
    setResult(event.target.value);
  };

  // console.log(active, delta, site, result, "[MyJobsControl]");
  return (
    <form style={{display:'block'}}
      className={classes.root} 
      noValidate 
      autoComplete="off">
      <div className={styles.FlexBox}>
        <div className={styles.TextFieldBox1}>  
          <TextField id={styles["text-box-days"]}
            size="small"
            label="DELTA-DAYS" 
            variant="outlined"
            value={delta}
            onChange={handleDelta}
            InputProps={{classes: {input: styles.Resize,},}}
            InputLabelProps={{classes: {root: styles.Resize,}}}
          />
        </div>
        <div className={styles.TextFieldBox2}>
          <TextField id={styles["text-box"]}
            size="small"
            select
            label="SITE"
            value={site}
            onChange={handleSite}
            variant="outlined"
            InputProps={{classes: {input: styles.Resize,},}}
            InputLabelProps={{classes: {root: styles.Resize,}}}
          >
          {sites.map((option) => (
            <MenuItem 
              className={styles.Resize} 
              key={option.value} 
              value={option.value}>
              {option.label} 
            </MenuItem>
          ))}
          </TextField>
        </div>
        <div className={styles.TextFieldBox3}>
          <TextField id={styles["text-box"]}
            size="small"
            select
            label="RESULT"
            value={result}
            onChange={handleResult}
            variant="outlined"
            InputProps={{classes: {input: styles.Resize,},}}
            InputLabelProps={{classes: {root: styles.Resize,}}}
          >
          {results.map((option) => (
            <MenuItem 
              className={styles.Resize} 
              key={option.value} 
              value={option.value}>
              {option.label}
            </MenuItem>
          ))}
          </TextField>
        </div>
        <div className={styles.SwitchBox}>
          <FormControlLabel
            className={styles.CheckBoxBlock}
            control={
              <Checkbox
                icon={<CheckBoxOutlineBlankIcon fontSize="small" />}
                checkedIcon={<CheckBoxIcon fontSize="small" />}
                checked={active}
                onChange={handleActive}
                name="active_only"
                color="primary"
                className={styles.CheckBox}
              />
            }
            label={<span className={styles.CheckBoxLabel}>
              active</span>}
          />
        </div>
        <div className={styles.ButtonBox}>
          <Button 
            size="large"
            variant="contained"
            color="primary"
            className={styles.Button}
            onClick={()=>{props.setAll(active, delta, site, result)}}
          >
            <span className={styles.CheckBoxLabel}>
              <b>SEARCH</b></span>
          </Button>
        </div>  
      </div>
    </form>
  );
}