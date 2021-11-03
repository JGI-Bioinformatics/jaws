import React, {Component } from 'react';

import ExpansionPanel from '@material-ui/core/ExpansionPanel';
import ExpansionPanelDetails from '@material-ui/core/ExpansionPanelDetails';
import ExpansionPanelSummary from '@material-ui/core/ExpansionPanelSummary';
import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import Divider from '@material-ui/core/Divider';

import Avatar from '@material-ui/core/Avatar';
import { Button, Tooltip, Box} from "@material-ui/core";
import Chip from '@material-ui/core/Chip';

// import {Redirect} from 'react-router-dom';
import {connect} from 'react-redux';


import Dialog from "../../UI/DraggableDialog/DraggableDialog";
import * as actions from "../../../store/actions";
import { avatars } from "../../../utils/avatars";

import Account from "./Account/Account";

// 
import classes from "./User.module.css"

class User extends Component {
  state = {
    dialog : false,
  }

  handleAvatarsDialog = (event) => {
    event.stopPropagation();
    // console.log("[handleAvatarsDialog]");
    this.setState({ dialog: true })
  }

  handleDialogClose = (idx) => {
    // console.log(idx, "[handleDialogClose");
    this.setState({dialog: false, });
    if (idx !== null){
      this.props.updateAvaterId(idx);
    }
  }

  render() {
    if (this.props.user === null){
      // return <Redirect to="/" />
      return <div style={{width: "100%", padding: "40px", textAlign: "center"}}>
        <h3>Not a logged user!</h3>
        </div>
    } 
    const useridentity = this.props.user.useridentity;

    //- User account summary banner
    const identities = useridentity.map( (itm, idx) => {
                              const avatar = <Avatar>{itm.organization.substr(0, 1)}</Avatar>

                              return <Chip
                                        key={idx}
                                        avatar={avatar}
                                        label={itm.organization}
                                        color="primary"
                                        variant="outlined"
                                      />
                        });

    let accountSammary = `${useridentity.length} Globus account`;
    accountSammary += useridentity.length > 1 ? "s:" : ":";

    // for avatar change dialog
    const avatarImgList = avatars.map( (img, idx) => (
          <Avatar key={idx} src={img}></Avatar> 
      ) );

    const dialog = this.state.dialog ? <Dialog 
                      click={this.handleDialogClose}
                      items={avatarImgList}
                      title="Change Avatar"
                      /> : null;
    
    return (
      <>
        <Divider variant="inset" component="hr" />
        <div className={classes.root}>
          <ExpansionPanel defaultExpanded={true}>
            <ExpansionPanelSummary className={classes.AccountSummaryExpansionPanel}
              expandIcon={<ExpandMoreIcon />}
              aria-controls="panel1c-content"
              id="panel1c-header"
            >
              <div className={classes.IdentityList}>
                <Tooltip title="Change Averta" arrow>
                  <Button onClick={ this.handleAvatarsDialog }>
                    <Avatar className={classes.LargeAvatar}
                            alt={this.props.user.userinfo.name} 
                            src={avatars[this.props.avatarId]}>
                    </Avatar>
                  </Button>
                </Tooltip>
                <Box className={classes.AccountSummary}>
                  {accountSammary} 
                  {identities}
                </Box>
              </div>
            </ExpansionPanelSummary>

            <ExpansionPanelDetails className={classes.details}>
              <Account />
            </ExpansionPanelDetails>
          </ExpansionPanel>
        </div>

        {dialog}
      </>
   );

  }

}

const mapStateToProps = state => {
  return {
    user: state.auth.user,
    avatarId: state.main.avatarId,
    userAccountLooked: state.userAccountLooked
  }
}

const mapDispatchToProps = dispatch => {
  return {
    updateAvaterId : (newId) => dispatch(actions.updateUserAvatarId(newId)),
  }
}

export default connect(mapStateToProps, mapDispatchToProps)(User);