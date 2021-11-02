import React, { useState, useRef } from "react";

import { PropTypes } from "prop-types";

import { withStyles } from "@material-ui/core/styles";
import Badge from "@material-ui/core/Badge";
import Avatar from "@material-ui/core/Avatar";
import NotificationsIcon from "@material-ui/icons/Notifications";
import ClickAwayListener from "@material-ui/core/ClickAwayListener";
import Grow from "@material-ui/core/Grow";
import Paper from "@material-ui/core/Paper";
import Popper from "@material-ui/core/Popper";
import MenuList from "@material-ui/core/MenuList";

import { avatars } from "../../utils/avatars";

import classes from "./PopMenu.module.css";

const badgeStyle = {
  right: 3,
  top: 10,
  padding: "0 4px",
  // marginRight: "20px",
};

const PendingBadge = withStyles((theme) => ({
  badge: {
    border: `2px solid ${theme.palette.background.paper}`,
    backgroundColor: "#f00",
    ...badgeStyle,
  },
}))(Badge);

const PopMenu = ({ content, callback, hasClosed, width, height }) => {
  const anchorRef = useRef();
  const [open, setOpen] = useState(false);

  const handleClose = (event) => {
    if (anchorRef.current && anchorRef.current.contains(event.target)) {
      return;
    }

    if (hasClosed) hasClosed();

    setOpen(false);
  };

  function handleListKeyDown(event) {
    if (event.key === "Tab") {
      event.preventDefault();
      setOpen(false);
    }
  }

  const getIcon = (data) => {
    let comp = null;

    switch (data.type) {
      case "avatar":
        comp = (
          <Avatar
            ref={anchorRef}
            className={data.css}
            onClick={() => {
              setOpen(!open);
            }}
            src={avatars[data.avatarId]}
          >
            {data.value} {/* src (icon) overwrite this feature, which is displaying the value string  */}
          </Avatar>
        );
        break;
      case "notice":
        comp = (
          <PendingBadge
            badgeContent={data.value}
            color="secondary"
            ref={anchorRef}
            className={classes.notificationBadge}
            onClick={() => {
              if (callback) callback();

              if (open) {
                if (hasClosed) hasClosed();
                setOpen(false);
              } else {
                setOpen(true);
              }
            }}
          >
            <NotificationsIcon className={classes.notificationAvatar} />
          </PendingBadge>
        );

        break;
      default:
        comp = null;
    }
    return comp;
  };

  const mainPanelStyle = {
    minWidth: `${width}px`,
    height: "auto",
    maxHeight: `${height}px`,
    overflow: "auto",
    zIndex: 20,
    border: "1px solid #d9d9d9",
    marginRight: "1.5em",
    outline: "none",
  };

  return (
    <>
      {getIcon(content.icon)}
      <Popper
        open={open}
        anchorEl={anchorRef.current}
        transition
        style={mainPanelStyle}
      >
        {({ TransitionProps, placement }) => (
          <Grow
            {...TransitionProps}
            style={{
              transformOrigin:
                placement === "bottom" ? "center top" : "center bottom",
            }}
          >
            <Paper>
              <ClickAwayListener onClickAway={handleClose}>
                <MenuList
                  autoFocusItem={open}
                  id="menu-list-grow"
                  onKeyDown={handleListKeyDown}
                >
                  {content.menus}
                </MenuList>
              </ClickAwayListener>
            </Paper>
          </Grow>
        )}
      </Popper>
    </>
  );
};

export default PopMenu;

PopMenu.defaultProps = {
  width: 300,
  height: 400,
  callback: null,
  hasClosed: null,
};

PopMenu.propTypes = {
  content: PropTypes.shape().isRequired,
  width: PropTypes.number,
  height: PropTypes.number,
  callback: PropTypes.func,
  hasClosed: PropTypes.func,
};
