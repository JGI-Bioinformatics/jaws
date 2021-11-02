import React, {Component} from 'react';

import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import Divider from '@material-ui/core/Divider';
import ListItemText from '@material-ui/core/ListItemText';
import ListItemAvatar from '@material-ui/core/ListItemAvatar';
import Avatar from '@material-ui/core/Avatar';
import Typography from '@material-ui/core/Typography';
import {Grid} from "@material-ui/core";

import {connect} from 'react-redux';

import classes from "./Account.module.css"

class Account extends Component {
  
  render() {
 
    const colorclass = [classes.orange, classes.purple, classes.blue, classes.gray];
    // const userId = this.props.userId;
    const useridentity = this.props.user.useridentity;

    // for multiple Globus account ...
    const carditems = useridentity
              .map( (itm, idx) => {
 
                  const title = [itm.name, itm.organization].join("  -  ")
                  return ( 
                    
                    <div key={itm.id}>
                      <Grid container spacing={2} justify="flex-start">
                        <Grid item md={6}  >
                          <Divider variant="inset" component="li" />
                          <ListItem alignItems="flex-start">
                            <ListItemAvatar>
                              <Avatar className={colorclass[idx%4]}>{itm.organization.substr(0, 1)}</Avatar>
                            </ListItemAvatar>
                            <ListItemText
                              primary={
                                title
                              }
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
                                  : {itm.email} <br />

                                  <Typography
                                    component="span"
                                    variant="body2"
                                    className={classes.inline}
                                    color="textPrimary"
                                  >
                                    username
                                  </Typography>
                                  : {itm.username} <br />

                                </React.Fragment>
                              }
                            />
                          </ListItem>
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
    user: state.auth.user,
  }
}

export default connect(mapStateToProps)(Account);