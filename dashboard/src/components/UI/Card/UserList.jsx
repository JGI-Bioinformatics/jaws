import React, {Component} from 'react';

import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import Divider from '@material-ui/core/Divider';
import ListItemText from '@material-ui/core/ListItemText';
import ListItemAvatar from '@material-ui/core/ListItemAvatar';
import Avatar from '@material-ui/core/Avatar';
import Typography from '@material-ui/core/Typography';
import {Grid} from "@material-ui/core";
// import Radio from '@material-ui/core/Radio';
// import RadioGroup from '@material-ui/core/RadioGroup';
// import FormControlLabel from '@material-ui/core/FormControlLabel';

import {Redirect} from 'react-router-dom';

import Checkbox from '@material-ui/core/Checkbox';
import FormControlLabel from '@material-ui/core/FormControlLabel';

import {connect} from 'react-redux';
import * as actionTypes from "../../../store/actions";

import classes from "./UserList.module.css"

class AlignItemsList extends Component {

  handleChange = (e, checked)=> {
    if (checked){
      const userId = e.target.value;
      this.props.updateUserId(userId);
    }
  }

  render() {
    if (this.props.user === null){
      return <Redirect to="/" />
    } 

    const colorclass = [classes.orange, classes.purple, classes.blue, classes.gray];
    const userId = this.props.userId;
    const useridentity = this.props.user.useridentity;

    const carditems = useridentity
              .map( (sec, idx) => {
                  const title = [sec.name, sec.organization].join("  -  ");
                  const avata = sec.organization.substr(0, 1);
                  
                  const cboxdefault = ( ( !userId && idx === 0 ) || (userId === sec.id) ) ? 
                        { checked:true} : { checked:false};
                  return ( 
                    
                    <div key={sec.id}>
                      <Grid container spacing={2} justify="flex-start">
                      <Grid item md={6}  >
                        <Divider variant="inset" component="li" />
                        <ListItem alignItems="flex-start">
                          <ListItemAvatar>
                            <Avatar className={colorclass[idx%4]}>{avata}</Avatar>
                          </ListItemAvatar>
                          <ListItemText
                            primary={title}
                            secondary={
                              <React.Fragment>
                                <Typography
                                  component="span"
                                  variant="body2"
                                  className={classes.inline}
                                  color="textPrimary"
                                >
                                  email
                                </Typography>
                                : {sec.email} <br />

                                <Typography
                                  component="span"
                                  variant="body2"
                                  className={classes.inline}
                                  color="textPrimary"
                                >
                                  username
                                </Typography>
                                : {sec.username} <br />
                                <Typography
                                  component="span"
                                  variant="body2"
                                  className={classes.inline}
                                  color="textPrimary"
                                >
                                  id
                                </Typography>
                                : {sec.id} <br />
                              </React.Fragment>
                            }
                          />
                        </ListItem>
                      </Grid>
                      <Grid item md={2}  >
                        
                        <FormControlLabel
                          control={<Checkbox
                            {...cboxdefault}
                            color="primary"
                            value={sec.id}
                            onChange={this.handleChange}
                            inputProps={{ 'aria-label': 'secondary checkbox' }}
                          /> }
                          label="in use"
                        />
                          
                        
                      </Grid>
                      </Grid>
                   
                    </div>
                      
                    )
                  }     
                )
  
    return (
      <List className={classes.root}>
          {carditems}
      </List>
    )
    
  }
}

const mapStateToProps = state => {
  return {
    userId : state.userId,
    user: state.user
  }
}

const mapDispatchToProps = dispatch => {
  return {
    updateUserId : (newId) => dispatch({ type:actionTypes.UPDATE_USER_ID, userId:newId})
  }
}

export default connect(mapStateToProps, mapDispatchToProps)(AlignItemsList);