import React, { Component } from "react";

// import InputLabel from "@material-ui/core/InputLabel";
import MenuItem from "@material-ui/core/MenuItem";
// import FormHelperText from "@material-ui/core/FormHelperText";
import FormControl from "@material-ui/core/FormControl";
import Select from "@material-ui/core/Select";
import { Grid, TextField, Button , InputLabel} from "@material-ui/core"

import classes from "./Search.module.css";


class Search extends Component {
  state = { 
      value: 0,
      search: "" 
  }

  searchTypeMap = [
    "Run ID",
    "WDL Name"
  ]

  handleSearch = () => {
    console.log(this.state, "[Search] handleSearch");
  }

  heandleSlectChange = (event) => {
    this.setState({ value: event.target.value});
    
  };

  handleSearchTextChange = (event) => {
    console.log(event.target.value, "[Search] handleSearchTextChange")
    this.setState({search: event.target.value})
  }

  render() { 
    const selections = this.searchTypeMap.map((e, idx) => <MenuItem key={e} value={idx}>{e}</MenuItem>)

    return ( 
    <div className="MainPanel">
      <h3> Site Search </h3>
      <Grid container spacing={1} className={classes.Master} justify="flex-start">
          <Grid item xs={12} sm={2} md={2}  className={classes.Item}>
            <FormControl id="search" className={classes.FormControl}>
              <InputLabel id="search-category">Search by</InputLabel>
              <Select
                labelId="search-category"
                id="demo-simple-select"
                value={this.state.value}
                onChange={this.heandleSlectChange}
              >
                {selections}
              </Select>
            </FormControl>
          </Grid>

          <Grid item xs={12} sm={6} md={7} className={classes.SearchVal}>
            <TextField label="Search value" 
              fullWidth 
              required
              defaultValue={this.state.search} 
              onChange={this.handleSearchTextChange}/>
          </Grid>

          <Grid item  xs={12} sm={4} md={2} className={classes.Button}>
            <Button 
            variant="contained" 
            size="small" 
            color="primary"
            onClick={this.handleSearch}>
              Do Search
            </Button>
          </Grid>
          
      </Grid>
    </div> );
  }
}
 
export default Search;