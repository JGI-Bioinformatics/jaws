import React, { useState} from "react";

import { PropTypes } from "prop-types";
import Avatar from "@material-ui/core/Avatar";
import {
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  Typography,
  Button, 
  Tooltip
} from "@material-ui/core";

import {connect} from 'react-redux';

import * as actions from "../../store/actions";
import { avatars } from "../../utils/avatars";
import Dialog from "../UI/DraggableDialog/DraggableDialog";

import classes from "./Header.module.css";

import PopMenu from "./PopMenu";

const UserAvatar = ({ 
  user, 
  handleLogout, 
  avatarId, 
  updateAvaterId, 
  handleUserAccountInfo
}) => {
  const [changeAvatarDialog, setChangeAvatarDialog] = useState(false);
 
  const handleChangeAvatar = ()=> {
    console.log("changeAvatar");
    setChangeAvatarDialog(true);
  }

  const handleChangeAvatarClose = (idx) => {
    // console.log(idx, "[handleDialogClose");
    setChangeAvatarDialog(false);
    if (idx !== null){
      updateAvaterId(idx);
    }
  }

  const fullname = `${user.first_name} ${user.last_name}`;
  const uname = fullname
    .split(" ")
    .map((e) => e[0])
    .join("")
    .toUpperCase();

  const content = {
    icon: {
      type: "avatar",
      value: uname,
      css: classes.userAvatar,
      avatarId: avatarId
    },
    menus: (
      <List>
        <ListItem divider>
          <ListItemIcon>
            <Tooltip title="Change Averta" arrow>
              <Button onClick={ handleChangeAvatar }>
                <Avatar className={classes.userAvatar}
                  alt={user.first_name} 
                  src={avatars[avatarId]}>
                </Avatar>
              </Button>
            </Tooltip>

            {/* <Avatar
              className={classes.userAvatar}
              alt={user.first_name} 
              src={avatars[avatarId]}
            /> */}
          </ListItemIcon>
          <ListItemText
            primary={
              <Typography className={classes.userMenuUserName}>
                {fullname}
              </Typography>
            }
            secondary={
              <>
                <Typography
                  component="span"
                  variant="body2"
                  className={classes.userMenuUserEmail}
                  style={{ color: "#444" }}
                >
                  {user.email_address.toUpperCase()}
                </Typography>
              </>
            }
          />
        </ListItem>
        <ListItem button onClick={handleUserAccountInfo}>
          <ListItemIcon> </ListItemIcon>
          <ListItemText
            primary={
              <Typography className={classes.userMenuText}>
                Account info
              </Typography>
            }
          />
        </ListItem>
        <ListItem button onClick={handleLogout}>
          <ListItemIcon> </ListItemIcon>
          <ListItemText
            primary={
              <Typography className={classes.userMenuText}>Logout</Typography>
            }
          />
        </ListItem>
      </List>
    ),
  };

  // for avatar change dialog
    const avatarImgList = avatars.map( (img, idx) => (
      <Avatar key={idx} src={img}></Avatar> 
  ) );

  const dialog = changeAvatarDialog ? <Dialog 
                  click={handleChangeAvatarClose}
                  items={avatarImgList}
                  title="Change Avatar"
                /> : null;

  return (
    <>
      <PopMenu content={content} width={200} />
      {dialog}
    </>
  );
};

const mapDispatchToProps = dispatch => {
  return {
    updateAvaterId : (newId) => dispatch(actions.updateUserAvatarId(newId)),
  }
}

export default connect(null, mapDispatchToProps)(UserAvatar);

UserAvatar.propTypes = {
  // handleLogout: PropTypes.func.isRequired,
  // honeycomb: PropTypes.shape({ sendUiInteractionSpan: PropTypes.func })
  //   .isRequired,
  avatarId: PropTypes.number.isRequired,
  handleLogout: PropTypes.func.isRequired,
  updateAvaterId: PropTypes.func.isRequired,
  user: PropTypes.shape({
    first_name: PropTypes.string.isRequired,
    last_name: PropTypes.string.isRequired,
    email_address: PropTypes.string.isRequired,
  }).isRequired,
};
